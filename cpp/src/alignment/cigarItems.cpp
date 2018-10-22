// All content Copyright (C) 2018 Genomics plc
#include "alignment/cigarItems.hpp"
#include <numeric>

#include "utils/logging.hpp"

#include "variant/type/variant.hpp"
#include "variant/snpFinder.hpp"

namespace wecall
{
namespace alignment
{

    referencePositions_t CigarMatch::getRelativeRefPositions( const int64_t & startRefPos ) const
    {
        referencePositions_t relativeRefPositions( int64_to_sizet( m_length ) );
        std::iota( relativeRefPositions.begin(), relativeRefPositions.end(), startRefPos );
        return relativeRefPositions;
    }

    //-----------------------------------------------------------------------------------------

    referencePositions_t CigarInsertion::getRelativeRefPositions( const int64_t & startRefPos ) const
    {
        return referencePositions_t( int64_to_sizet( m_length ), emptyPos );
    }

    //-----------------------------------------------------------------------------------------

    referencePositions_t CigarSoftClip::getRelativeRefPositions( const int64_t & startRefPos ) const
    {
        return referencePositions_t( int64_to_sizet( m_length ), emptyPos );
    }

    //-----------------------------------------------------------------------------------------

    void CigarItem::moveOffsets( offsetsPtr_t offsets ) const
    {
        offsets->read += this->lengthInSeq();
        offsets->ref += this->lengthInRef();
    }

    //-----------------------------------------------------------------------------------------

    variant::variantSet_t CigarMatch::getVariants( variant::variantGenerationDataPtr_t variantGenerationData,
                                                   offsetsPtr_t offsets ) const
    {
        if ( ( variantGenerationData->readStartPos + offsets->ref + m_length ) <
             variantGenerationData->refSeq->start() )
        {
            return {};
        }
        // else
        return variant::SNPFinder( variantGenerationData ).findSNPsInReadSegment( offsets, m_length );
    }

    //-----------------------------------------------------------------------------------------

    variant::variantSet_t CigarInsertion::getVariants( variant::variantGenerationDataPtr_t varGenData,
                                                       offsetsPtr_t offsets ) const
    {
        if ( ( varGenData->readStartPos + offsets->ref ) < varGenData->refSeq->start() )
        {
            return {};
        }

        if ( offsets->read <= 0 or
             static_cast< std::size_t >( offsets->read + m_length ) >= varGenData->readSeq.size() )
        {
            return {};
        }

        if ( m_length == 0 )
        {
            return {};
        }

        const auto added = varGenData->readSeq.substr( offsets->read, m_length );

        // Note: readStart + refOffset - 1 below because reference position is the base before the insertion
        caller::Region region( varGenData->refSeq->contig(), varGenData->readStartPos + offsets->ref,
                               varGenData->readStartPos + offsets->ref );

        variant::varPtr_t insertion = std::make_shared< variant::Variant >( varGenData->refSeq, region, added, false );

        return {insertion};
    }

    //-----------------------------------------------------------------------------------------

    variant::variantSet_t CigarDeletion::getVariants( variant::variantGenerationDataPtr_t varGenData,
                                                      offsetsPtr_t offsets ) const
    {
        if ( ( varGenData->readStartPos + offsets->ref ) < varGenData->refSeq->start() )
        {
            return {};
        }

        if ( offsets->read <= 0 or static_cast< std::size_t >( offsets->read ) >= varGenData->readSeq.size() )
        {
            return {};
        }

        // else
        const caller::Region deletionLocation( varGenData->refSeq->contig(), varGenData->readStartPos + offsets->ref,
                                               varGenData->readStartPos + offsets->ref + m_length );
        // Deletion position here is last ref base before deletion
        variant::varPtr_t deletion =
            std::make_shared< variant::Variant >( varGenData->refSeq, deletionLocation, "", false );

        return {deletion};
    }

}  // namespace variant
}  // namespace wecall
