// All content Copyright (C) 2018 Genomics plc
#include <iomanip>
#include <boost/algorithm/string/join.hpp>

#include "vcf/record.hpp"
#include "variant/type/variant.hpp"
#include "utils/write.hpp"

namespace echidna
{
namespace vcf
{
    Record::Record( std::string contig,
                    std::size_t pos,
                    std::set< std::string > ids,
                    std::string ref,
                    std::vector< std::string > alts,
                    phred_t qual,
                    std::set< std::string > filters,
                    Info info,
                    SampleInfo sampleInfo )
        : m_contig( contig ),
          m_pos( pos ),
          m_ids( ids ),
          m_ref( ref ),
          m_alts( alts ),
          m_qual( qual ),
          m_filters( filters ),
          m_info( info ),
          m_sampleInfo( sampleInfo )
    {
        // Nothing to do here
    }

    std::vector< variant::varPtr_t > Record::getVariants(
        const utils::referenceSequencePtr_t & referenceSequence ) const
    {
        std::vector< variant::varPtr_t > variants;
        for ( auto alt : m_alts )
        {
            auto internalPos = m_pos - 1;

            const caller::Region refRegion( m_contig, internalPos, internalPos + m_ref.size() );
            const auto var = std::make_shared< variant::Variant >( referenceSequence, refRegion, alt, false );
            if ( var->refSequence().sequence() == m_ref )
            {
                const auto trimmed = var->getTrimmed();
                if ( var->sequence() != trimmed->sequence() and m_ref.size() > 1 and alt.size() > 1 )
                {
                    ECHIDNA_LOG( WARNING, "Trimming input variant from: " + var->toString() + " to " +
                                              trimmed->toString() + "." );
                }

                variants.push_back( trimmed );
            }
            else
            {
                ECHIDNA_LOG( WARNING, "Ignoring record at: " << var->region() << "due to incompatible reference." );
            }
        }
        return variants;
    }

    std::ostream & operator<<( std::ostream & out, const Record & record )
    {
        const auto sampleInfoItemSeparator = ":";

        out << record.m_contig << constants::vcfRecordColumnSeparator;
        out << record.m_pos << constants::vcfRecordColumnSeparator;
        out << ( record.m_ids.empty() ? constants::vcfUnknownValue : boost::algorithm::join( record.m_ids, ";" ) )
            << constants::vcfRecordColumnSeparator;
        out << record.m_ref << constants::vcfRecordColumnSeparator;
        out << boost::algorithm::join( record.m_alts, "," ) << constants::vcfRecordColumnSeparator;
        out << caller::serialise_phred( record.m_qual ) << constants::vcfRecordColumnSeparator;
        out << ( record.m_filters.empty() ? "PASS" : boost::algorithm::join( record.m_filters, ";" ) )
            << constants::vcfRecordColumnSeparator;

        // Output info ID, value pairs
        std::vector< std::string > infoItems;
        for ( auto & idValuePair : record.m_info )
        {
            infoItems.push_back( idValuePair.first + "=" + boost::algorithm::join( idValuePair.second, "," ) );
        }
        out << boost::algorithm::join( infoItems, ";" ) << constants::vcfRecordColumnSeparator;

        // Output sample info IDs first...
        std::vector< std::string > sampleColumns = {
            boost::algorithm::join( record.m_sampleInfo.first, sampleInfoItemSeparator )};

        // ...then a matching set of values for each sample (tab separated)
        for ( auto & sampleInfoValues : record.m_sampleInfo.second )
        {
            std::vector< std::string > innerValues;
            for ( auto it = sampleInfoValues.begin(); it != sampleInfoValues.end(); ++it )
            {
                innerValues.push_back( boost::algorithm::join( *it, "," ) );
            }

            sampleColumns.push_back( boost::algorithm::join( innerValues, sampleInfoItemSeparator ) );
        }
        out << boost::algorithm::join( sampleColumns, constants::vcfRecordColumnSeparator ) << std::endl;
        return out;
    }
}
}
