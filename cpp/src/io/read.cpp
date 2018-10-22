// All content Copyright (C) 2018 Genomics plc
#include "io/read.hpp"
#include "utils/exceptions.hpp"
#include "utils/logging.hpp"
#include "alignment/cigar.hpp"
#include "alignment/cigarItems.hpp"
#include "io/pysam.hpp"
#include "read.hpp"

#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <utility>
#include <cassert>
#include <algorithm>
#include <sstream>
#include <utils/sequence.hpp>
#include <samtools/bam.h>

namespace wecall
{
namespace io
{
    using alignment::cigarFlags;

    uint8_t * pysam_bam1_qname( const bam1_t * b ) { return (uint8_t *)( b->data ); }
    // Look-up table for BAM sequence characters
    static const char * bamTable = "=ACMGRSVTWYHKDBN";

    //-----------------------------------------------------------------------------------------

    Read::ReadParams::ReadParams( const bam1_t * bamRecord )
        : m_qualities( int32_to_sizet( bamRecord->core.l_qseq ), 0 ),
          m_cigar( bamRecord->core.n_cigar, pysam_bam1_cigar( bamRecord ) ),
          m_tid( bamRecord->core.tid ),
          m_startPos( bamRecord->core.pos ),
          m_endPos( m_startPos + bamRecord->core.l_qseq ),
          m_flag( bamRecord->core.flag ),
          m_mappingQuality( bamRecord->core.qual ),
          m_insertSize( bamRecord->core.isize ),
          m_mateTid( bamRecord->core.mtid ),
          m_mateStartPos( bamRecord->core.mpos )
    {

        const unsigned char * qname_start = pysam_bam1_qname( bamRecord );
        m_qname.assign( qname_start, qname_start + bamRecord->core.l_qname - 1 );

        const uint8_t * readGroupTag = bam_aux_get( bamRecord, "RG" );

        if ( readGroupTag )
        {
            m_readGroupID = bam_aux2Z( readGroupTag );
        }

        const unsigned char * p = pysam_bam1_seq( bamRecord );
        const unsigned char * q = pysam_bam1_qual( bamRecord );

        std::string sequenceStr( int32_to_sizet( bamRecord->core.l_qseq ), ' ' );

        auto end = sequenceStr.end();
        auto evenEnd = sequenceStr.size() % 2 ? end - 1 : end;
        auto seqIt = sequenceStr.begin();
        std::copy( q, q + sequenceStr.size(), m_qualities.begin() );

        for ( ; seqIt != evenEnd; p++ )
        {
            *( seqIt++ ) = bamTable[( *p >> 4 ) & 0xf];
            *( seqIt++ ) = bamTable[*p & 0xf];
        }
        if ( seqIt != end )
            *seqIt = bamTable[( *p >> 4 ) & 0xf];

        m_sequenceStr = sequenceStr;
        if ( sequenceStr.size() == 0 || q[0] == 0xff )
        {
            throw utils::wecall_exception( "Invalid bamRecord in Read constructor" );
        }
    }

    Read::ReadParams::ReadParams( std::string seq,
                                  utils::QualitySequence qual,
                                  std::string readGroupID,
                                  alignment::Cigar cigar,
                                  int32_t tid,
                                  int32_t startPos,
                                  uint16_t flag,
                                  uint8_t mappingQuality,
                                  int32_t insertSize,
                                  int32_t mateTid,
                                  int32_t mateStartPos,
                                  const std::string & qname )
        : m_sequenceStr( seq ),
          m_qualities( qual ),
          m_qname( qname ),
          m_readGroupID( readGroupID ),
          m_cigar( cigar ),
          m_tid( tid ),
          m_startPos( startPos ),
          m_endPos( startPos + seq.size() ),
          m_flag( flag ),
          m_mappingQuality( mappingQuality ),
          m_insertSize( insertSize ),
          m_mateTid( mateTid ),
          m_mateStartPos( mateStartPos )
    {
    }

    Read::ReadParams::ReadParams()
        : Read::ReadParams( std::string(),
                            utils::QualitySequence( "" ),
                            std::string(),
                            alignment::Cigar( "" ),
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            std::string() )
    {
    }

    Read::Read( const Read::ReadParams & params, utils::referenceSequencePtr_t refSequence )
        : m_sequence(),
          m_qualities( params.m_qualities ),
          m_qname( params.m_qname ),
          m_readGroupID( params.m_readGroupID ),
          m_cigar( params.m_cigar ),
          m_tid( params.m_tid ),
          m_startPos( params.m_startPos ),
          m_endPos( params.m_endPos ),
          m_alignedEndPos( m_startPos + m_cigar.lengthInRef() ),
          m_insertSize( params.m_insertSize ),
          m_mateTid( params.m_mateTid ),
          m_mateStartPos( params.m_mateStartPos ),
          m_refSequence( refSequence ),
          m_flag( params.m_flag ),
          m_mappingQuality( params.m_mappingQuality )
    {
        if ( m_endPos == m_alignedEndPos )
        {
            auto ref_sequence = getRefSequenceRange();
            m_isReference =
                std::equal( params.m_sequenceStr.cbegin(), params.m_sequenceStr.cend(), ref_sequence.first );
        }

        if ( not m_isReference )
        {
            const_cast< utils::BasePairSequence & >( m_sequence ) = params.m_sequenceStr;
        }
    }

    Read::Read( const bam1_t * bamRecord, utils::referenceSequencePtr_t refSequence )
        : Read( Read::ReadParams( bamRecord ), refSequence )
    {
    }

    utils::BasePairSequence Read::makeRefSequence() const
    {
        return m_refSequence->subseq( caller::Region( m_refSequence->contig(), m_startPos, m_alignedEndPos ) )
            .sequence();
    }

    std::pair< utils::BasePairSequence::const_iterator, utils::BasePairSequence::const_iterator >
    Read::getRefSequenceRange() const
    {
        return m_refSequence->subseqRange( caller::Region( m_refSequence->contig(), m_startPos, m_alignedEndPos ) );
    }

    //-----------------------------------------------------------------------------------------

    Read::Read( utils::BasePairSequence seq,
                utils::QualitySequence qual,
                std::string readGroupID,
                alignment::Cigar cigar,
                int32_t tid,
                int32_t startPos,
                uint16_t flag,
                uint8_t mappingQuality,
                int32_t insertSize,
                int32_t mateTid,
                int32_t mateStartPos,
                utils::referenceSequencePtr_t refSequence,
                std::string qname )
        : Read( Read::ReadParams( seq.str(),
                                  qual,
                                  readGroupID,
                                  cigar,
                                  tid,
                                  startPos,
                                  flag,
                                  mappingQuality,
                                  insertSize,
                                  mateTid,
                                  mateStartPos,
                                  qname ),
                refSequence )
    {
        ECHIDNA_ASSERT( seq.size() == m_qualities.size(), "Read constructed with inconsistent quality string " +
                                                              std::to_string( m_sequence.size() ) + " != " +
                                                              std::to_string( m_qualities.size() ) );
        ECHIDNA_ASSERT( seq.size() == int64_to_sizet( cigar.lengthInSeq() ),
                        "Read constructed with inconsistent cigar " + std::to_string( m_sequence.size() ) + " != " +
                            std::to_string( cigar.lengthInSeq() ) + " " + cigar.toString() );
    }

    //-----------------------------------------------------------------------------------------

    std::string Read::toString() const
    {
        std::stringstream repr;

        repr << m_startPos << " --> " << m_alignedEndPos << "  cigar= " << m_cigar.toString() << " " << getQName()
             << " reference=" << std::to_string( m_isReference ) << std::endl;
        repr << this->sequence() << "  insert size = " << m_insertSize << std::endl;
        for ( std::size_t i = 0; i < m_qualities.size(); ++i )
        {
            repr << static_cast< char >( m_qualities[i] + 33 );
        }
        repr << "  mapQ = " << std::to_string( m_mappingQuality );

        return repr.str();
    }

    utils::Interval Read::getMaximalReadInterval() const
    {
        const auto readRefStartPos = this->getStartPos();
        const auto readRefEndPos = this->getAlignedEndPos();

        const auto lengthBeforeAlignedStartPos = this->getLengthBeforeAlignedStartPos();
        const auto lengthAfterAlignedEndPos = this->getLengthAfterAlignedEndPos();

        const auto adjustedStartPos = readRefStartPos - lengthBeforeAlignedStartPos;
        const auto adjustedEndPos = readRefEndPos + lengthAfterAlignedEndPos;

        return utils::Interval( adjustedStartPos, adjustedEndPos );
    }

    //-----------------------------------------------------------------------------------------

    int64_t Read::getLengthAfterAlignedEndPos() const { return m_cigar.lengthAfterRefEndPos(); }

    //-----------------------------------------------------------------------------------------

    bool Read::isReverse() const { return ( ( m_flag & BAM_FREVERSE ) != 0 ); }

    //-----------------------------------------------------------------------------------------

    bool Read::isPaired() const { return ( ( m_flag & BAM_FPAIRED ) != 0 ); }

    //-----------------------------------------------------------------------------------------

    bool Read::isDuplicate() const { return ( ( m_flag & BAM_FDUP ) != 0 ); }

    //-----------------------------------------------------------------------------------------

    bool Read::isSecondary() const { return ( ( m_flag & BAM_FSECONDARY ) != 0 ); }

    //-----------------------------------------------------------------------------------------

    bool Read::isUnMapped() const { return ( ( m_flag & BAM_FUNMAP ) != 0 ); }

    //-----------------------------------------------------------------------------------------

    bool Read::isMateUnMapped() const { return ( ( m_flag & BAM_FMUNMAP ) != 0 ); }

    //-----------------------------------------------------------------------------------------

    bool Read::isMateReverse() const { return ( ( m_flag & BAM_FMREVERSE ) != 0 ); }

    //-----------------------------------------------------------------------------------------

    bool Read::isReadOne() const { return ( ( m_flag & BAM_FREAD1 ) != 0 ); }

    //-----------------------------------------------------------------------------------------

    bool Read::isProperPair() const { return ( ( m_flag & BAM_FPROPER_PAIR ) != 0 ); }

    //-----------------------------------------------------------------------------------------

    void Read::trimOverlap()
    {
        const auto insertSize = getInsertSize();
        if ( isReadOne() or not isProperPair() or insertSize == 0 )
        {
            return;
        }
        // else
        const auto absIns = std::abs( insertSize );
        // ES: Assume that the reads have the same length. Ideally not do this.
        const auto overlapLength = 2 * getLength() - absIns;

        if ( isReverse() )
        {
            for ( auto qualIt = m_qualities.begin();
                  qualIt != m_qualities.end() and qualIt < m_qualities.begin() + overlapLength; ++qualIt )
            {
                *qualIt = constants::minAllowedQualityScore;
            }
        }
        else
        {
            for ( auto qualIt = m_qualities.rbegin();
                  qualIt != m_qualities.rend() and qualIt < m_qualities.rbegin() + overlapLength; ++qualIt )
            {
                *qualIt = constants::minAllowedQualityScore;
            }
        }
    }

    //-----------------------------------------------------------------------------------------

    void Read::trimReadOfShortFragment()
    {
        const auto insertSize = getInsertSize();

        if ( not isProperPair() or insertSize == 0 )
        {
            return;
        }

        // else
        const auto absIns = std::abs( insertSize );

        if ( absIns > this->getLength() )
        {
            return;
        }

        if ( isReverse() )
        {
            for ( auto qualIt = m_qualities.begin(); qualIt < m_qualities.end() - absIns; ++qualIt )
            {
                *qualIt = constants::minAllowedQualityScore;
            }
        }

        else
        {
            for ( auto qualIt = m_qualities.rbegin(); qualIt < m_qualities.rend() - absIns; ++qualIt )
            {
                *qualIt = constants::minAllowedQualityScore;
            }
        }
    }

    //-----------------------------------------------------------------------------------------

    alignment::referencePositions_t Read::getReferencePositions() const
    {
        return m_cigar.getRefPositions( m_startPos );
    }

    utils::Interval Read::getIntervalInRead( const utils::Interval & refInterval ) const
    {
        utils::Interval relativeInterval = refInterval - this->getStartPos();
        if ( relativeInterval.end() <= 0L )
        {
            relativeInterval = utils::Interval( 0L, 0L );
        }
        else if ( relativeInterval.start() >= this->getAlignedLength() )
        {
            relativeInterval = utils::Interval( this->getAlignedLength(), this->getAlignedLength() );
        }
        else
        {
            relativeInterval = relativeInterval.getIntersect( utils::Interval( 0L, this->getAlignedLength() ) );
        }

        return cigar().getInverseInterval( relativeInterval );
    }

    std::vector< variant::varPtr_t > Read::getVariants() const
    {
        if ( this->isReference() )
        {
            return {};
        }

        const auto variantGenerationData = std::make_shared< const variant::VariantGenerationData >(
            m_refSequence, this->getStartPos(), this->sequence() );

        auto offsets = std::make_shared< alignment::Offsets >( 0L, 0L );

        std::vector< variant::varPtr_t > variants;

        for ( auto cigarIt = m_cigar.cbegin(); cigarIt != m_cigar.cend(); ++cigarIt )
        {
            const auto cigarVariants = ( *cigarIt )->getVariants( variantGenerationData, offsets );
            variants.insert( variants.end(), cigarVariants.begin(), cigarVariants.end() );
            ( *cigarIt )->moveOffsets( offsets );
        }

        return variants;
    }

    std::vector< variant::breakpointPtr_t > Read::getBreakpoints() const
    {
        if ( m_cigar.size() < 2 )
        {
            return {};
        }

        const auto contig = m_refSequence->contig();
        std::vector< variant::breakpointPtr_t > breakpoints;
        if ( m_cigar.front()->isSoftClipped() )
        {
            utils::BasePairSequence sequence = this->sequence().substr( 0, m_cigar.front()->lengthInSeq() );
            breakpoints.push_back( std::make_shared< variant::Breakpoint >( contig, m_startPos, false, sequence ) );
        }
        else if ( m_cigar.front()->isHardClipped() )
        {
            breakpoints.push_back( std::make_shared< variant::Breakpoint >( contig, m_startPos, false, "" ) );
        }

        if ( m_cigar.back()->isSoftClipped() )
        {
            utils::BasePairSequence sequence = this->sequence().substr(
                m_qualities.size() - m_cigar.back()->lengthInSeq(), m_cigar.back()->lengthInSeq() );
            breakpoints.push_back( std::make_shared< variant::Breakpoint >( contig, m_alignedEndPos, true, sequence ) );
        }
        else if ( m_cigar.back()->isHardClipped() )
        {
            breakpoints.push_back( std::make_shared< variant::Breakpoint >( contig, m_alignedEndPos, true, "" ) );
        }

        for ( auto breakpoint : breakpoints )
        {
            if ( not this->isMateUnMapped() and this->isMateOnSameContig() )
            {
                breakpoint->addQueryRegion( caller::Region( contig, this->getMateIntervalInRef() ) );
            }
        }

        return breakpoints;
    }
}
}
