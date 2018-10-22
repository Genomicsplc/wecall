// All content Copyright (C) 2018 Genomics plc
#include "variant/snpFinder.hpp"
#include "variant/type/variant.hpp"
#include "variant/clustering.hpp"

namespace wecall
{
namespace variant
{

    int64_t SNPFinder::refIndexFromReadIndex( const int64_t readIndex, const alignment::offsetsPtr_t offsets ) const
    {
        return readIndex - offsets->read + offsets->ref + m_varGenData->readStartPos;
    }

    int64_t SNPFinder::readIndexFromRefIndexInStr( const alignment::offsetsPtr_t offsets ) const
    {
        return offsets->read - offsets->ref - m_varGenData->readStartPos + m_varGenData->refSeq->start();
    }

    variantSet_t SNPFinder::findSNPsInReadSegment( const alignment::offsetsPtr_t offsets, const int64_t length ) const
    {
        variantSet_t snps;

        const int64_t endReadIndex =
            std::min( offsets->read + length, static_cast< int64_t >( m_varGenData->readSeq.size() ) );

        int64_t readIndex = std::max( offsets->read, readIndexFromRefIndexInStr( offsets ) );
        int64_t refIndex = refIndexFromReadIndex( readIndex, offsets );

        if ( m_varGenData->refSeq->start() > refIndex )
        {
            const auto diffFromRegionStart = m_varGenData->refSeq->start() - refIndex;
            readIndex += diffFromRegionStart;
            refIndex += diffFromRegionStart;
        }

        for ( ; readIndex < endReadIndex and refIndex < m_varGenData->refSeq->end(); ++readIndex, ++refIndex )
        {
            const auto refChar = m_varGenData->refSeq->at( refIndex );
            const auto readChar = m_varGenData->readSeq[readIndex];

            if ( refChar != m_varGenData->readSeq[readIndex] )
            {
                const caller::Region snpRegion( m_varGenData->refSeq->contig(), refIndex, refIndex + 1 );
                snps.insert( std::make_shared< Variant >( m_varGenData->refSeq, snpRegion,
                                                          utils::BasePairSequence( 1, readChar ), true ) );
            }
        }
        return snps;
    }
}
}
