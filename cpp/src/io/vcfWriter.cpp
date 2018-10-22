// All content Copyright (C) 2018 Genomics plc
#include "io/vcfWriter.hpp"

#include <boost/filesystem.hpp>
#include "utils/timer.hpp"

#include <stdexcept>
#include <algorithm>

namespace wecall
{
using variant::genotypePtr_t;
using variant::varPtr_t;
using caller::Annotation;

namespace io
{
    //-----------------------------------------------------------------------------------------

    VCFWriter::VCFWriter( const std::string & outputFilename, bool outputRefCalls, bool outputPhasedGenotypes )
        : m_headerWritten( false ),
          m_outputRefCalls( outputRefCalls ),
          m_outputPhasedGenotypes( outputPhasedGenotypes ),
          m_file( outputFilename.c_str() ),
          m_contig( "" ),
          m_timer( std::make_shared< utils::Timer >( "IO", utils::fileMetaData( outputFilename ) ) )
    {
        if ( not( m_file.is_open() ) )
        {
            throw utils::wecall_exception( "Could not open VCF file for writing" );
        }
    }

    //-----------------------------------------------------------------------------------------

    VCFWriter::~VCFWriter() { m_file.close(); }

    //-----------------------------------------------------------------------------------------

    void VCFWriter::writeHeader( const std::string & userSpecifiedFormat,
                                 const caller::params::Application & applicationParams,
                                 const std::string & ref,
                                 const std::vector< std::string > & sampleNames,
                                 std::vector< vcf::FilterDesc > filterDescs,
                                 std::vector< caller::Region > contigs )
    {
        utils::ScopedTimerTrigger scopedTimerTrigger( m_timer );
        if ( m_headerWritten )
        {
            throw utils::wecall_exception( "Tried to write VCF header twice" );
        }

        // TODO - This should be tied in with the list in Annotation.cpp/hpp

        const std::string shortRefName = boost::filesystem::path( ref ).stem().string();

        const vcf::Header header( userSpecifiedFormat, applicationParams, shortRefName, sampleNames,
                                  vcf::info::getVCFKeys( m_outputRefCalls ),
                                  vcf::format::getVCFKeys( m_outputPhasedGenotypes, m_outputRefCalls ), filterDescs,
                                  contigs );

        m_file << header;
        m_headerWritten = true;
    }

    //-----------------------------------------------------------------------------------------

    void VCFWriter::contig( const std::string & contig ) { m_contig = contig; }

    //-----------------------------------------------------------------------------------------

    bool hasCanonicalBases( const std::string & ref, const std::string & alt )
    {
        const auto baseNonCanonical = []( char c )
        {
            return c != 'A' and c != 'C' and c != 'G' and c != 'T';
        };

        if ( alt != constants::vcfRefAltValue )
        {
            if ( std::find_if( alt.cbegin(), alt.cend(), baseNonCanonical ) != alt.cend() )
            {
                return false;
            }
        }
        return std::find_if( ref.cbegin(), ref.cend(), baseNonCanonical ) == ref.cend();
    }

    void VCFWriter::writeCallSet( io::FastaFile & refFile, const caller::callVector_t & calls )
    {
        utils::ScopedTimerTrigger scopedTimerTrigger( m_timer );
        if ( not m_headerWritten )
        {
            throw utils::wecall_exception( "Tried to write VCF record before header" );
        }

        for ( callIt_t it = calls.cbegin(); it != calls.cend(); ++it )
        {

            const auto pos = it->interval.start();
            std::string ref;
            std::string alt;

            if ( it->isRefCall() )
            {
                // weCall specific VCF notation for a reference call
                ref = refFile.getSequence( caller::Region( m_contig, pos, pos + 1 ) ).sequence().str();
                alt = constants::vcfRefAltValue;
            }
            else
            {
                if ( it->var->contig() != m_contig )
                {
                    std::ostringstream message;
                    message << "contigs don't match (expected: '" << m_contig << "'; got '" << it->var->contig()
                            << "') for variant at " << it->var->toString();
                    throw utils::wecall_exception( message.str().c_str() );
                }
                const auto refAndAlts = compileRefsAndAlts( it->var, refFile );
                ref = refAndAlts.first.sequence().str();
                alt = refAndAlts.second;
            }

            if ( hasCanonicalBases( ref, alt ) )
            {
                vcf::Record record( m_contig, pos + 1, it->varIds, ref, {alt}, it->qual, it->filters, compileInfo( it ),
                                    compileSampleInfo( it, m_outputPhasedGenotypes ) );
                m_file << record;
            }
            else
            {
                ECHIDNA_LOG( DEBUG, "Not outputting due to non-canonical bases:\t" << m_contig << "\t" << pos + 1
                                                                                   << "\t" << ref << "\t" << alt );
            }
        }
    }

    //-----------------------------------------------------------------------------------------

    const std::pair< utils::ReferenceSequence, std::string > VCFWriter::compileRefsAndAlts(
        variant::varPtr_t var,
        const io::FastaFile & refFile )
    {
        // Collect VCF formatted refs and alts
        const auto contig = var->contig();
        const int64_t minPos = var->zeroIndexedVcfPosition();
        const auto refSequence = refFile.getSequence( caller::Region( contig, minPos, var->end() ) );

        const caller::Region leftRegion( contig, minPos, var->start() );
        const auto paddedAlt = refSequence.subseq( leftRegion ).sequence() + var->sequence();

        return make_pair( refSequence, paddedAlt.str() );
    }

    //-----------------------------------------------------------------------------------------

    const vcf::Info VCFWriter::compileInfo( callIt_t it )
    {
        vcf::Info info;
        for ( AnnotationPtr_t annotation : it->getAnnotations() )
        {
            info.emplace_back( std::make_pair( annotation->getID(), annotation->getValues() ) );
        }
        return info;
    }

    //-----------------------------------------------------------------------------------------

    const vcf::SampleInfo VCFWriter::compileSampleInfo( callIt_t it, const bool outputPhasedGenotypes )
    {
        // Start with the genotype to variant mappings

        vcf::SampleInfoFormat format = {vcf::format::GT_key};
        // Extract format only for first variant & sample (must be same for all)
        for ( const auto & annotation : it->samples[0].getAnnotations() )
        {
            format.emplace_back( annotation->getID() );
        }

        // Assume same annotations for all variants & samples
        std::vector< vcf::SampleInfoValues > values = compileGenVarMap( it, outputPhasedGenotypes );

        for ( std::size_t sampleIndex = 0; sampleIndex < it->samples.size(); ++sampleIndex )
        {
            auto & sampleInfoValue = values[sampleIndex];
            for ( const auto & annotation : it->samples[sampleIndex].getAnnotations() )
            {
                sampleInfoValue.emplace_back( annotation->getValues() );
            }
        }

        return std::make_pair( format, values );
    }

    //-----------------------------------------------------------------------------------------

    const std::vector< vcf::SampleInfoValues > VCFWriter::compileGenVarMap( callIt_t it,
                                                                            const bool outputPhasedGenotypes )
    {
        constexpr int UNKNOWN_CALL = -1;
        constexpr int VARIANT_CALL = +1;

        // Initialise genotype to variant map with zeros (=> REF)

        const auto nSamples = it->samples.size();
        std::vector< std::vector< int > > genVarMap( nSamples );

        for ( std::size_t s = 0; s < nSamples; ++s )
        {
            genVarMap[s] = std::vector< int >( it->samples[s].genotypeCalls.size(), 0 );
        }

        // Replace with ALT index (starting with 1) where variant called or
        // the VCF unknown value where not known.

        for ( std::size_t s = 0; s < nSamples; ++s )
        {
            const genoCall_t & called = it->samples[s].genotypeCalls;
            auto & genotypeThisSample = genVarMap[s];
            for ( std::size_t i = 0; i < genotypeThisSample.size(); ++i )
            {
                switch ( called[i] )
                {
                case caller::Call::VAR:
                    genotypeThisSample[i] = VARIANT_CALL;
                    break;
                case caller::Call::UNKNOWN:
                    genotypeThisSample[i] = UNKNOWN_CALL;
                    break;
                case caller::Call::REF:
                    break;
                }
            }
            if ( genotypeThisSample.empty() )
            {
                genotypeThisSample = {UNKNOWN_CALL};
            }
        }

        std::vector< vcf::SampleInfoValues > values( nSamples );

        for ( std::size_t sampleIndex = 0; sampleIndex < nSamples; ++sampleIndex )
        {
            auto genotypeValues = genVarMap[sampleIndex];
            if ( not outputPhasedGenotypes )
            {
                std::sort( genotypeValues.begin(), genotypeValues.end() );
            }

            const auto serialise_integral_type = []( const int val )
            {
                return val == UNKNOWN_CALL ? constants::vcfUnknownValue : std::to_string( val );
            };

            std::vector< std::string > serialiasedGenotypes( genotypeValues.size() );
            std::transform( genotypeValues.begin(), genotypeValues.end(), serialiasedGenotypes.begin(),
                            serialise_integral_type );

            const auto genotypeSeparator = outputPhasedGenotypes ? constants::vcfPhasedGenotypeDeliminator
                                                                 : constants::vcfUnphasedGenotypeDeliminator;
            const std::string gtValue = boost::algorithm::join( serialiasedGenotypes, genotypeSeparator );
            values[sampleIndex] = {{gtValue}};
        }

        return values;
    }

    //-----------------------------------------------------------------------------------------
}
}
