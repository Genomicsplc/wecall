// All content Copyright (C) 2018 Genomics plc
#include "variant/variantNormalizer.hpp"
//#include <boost/numeric/ublas/io.hpp>
#include <stack>
#include "utils/interval.hpp"

namespace wecall
{
namespace variant
{
    variantSet_t VariantNormalizer::getNormalized( const caller::Region & paddedRegion,
                                                   const utils::BasePairSequence & sequence,
                                                   const boost::optional< variantSet_t > unnormalizedVariants ) const
    {
        caller::Region region = paddedRegion;
        utils::BasePairSequence referenceString = m_referenceSequence->subseq( paddedRegion ).sequence();
        utils::BasePairSequence altString = sequence;

        if ( unnormalizedVariants.is_initialized() )
        {
            const auto variants = unnormalizedVariants.get();

            if ( variants.empty() )
            {
                return {};
            }

            const auto startComp = []( const varPtr_t & a, const varPtr_t & b )
            {
                return a->start() < b->start();
            };
            const auto minElement = std::min_element( variants.cbegin(), variants.cend(), startComp );
            const auto start = ( *minElement )->start();

            const auto endComp = []( const varPtr_t & a, const varPtr_t & b )
            {
                return a->end() < b->end();
            };
            const auto maxElement = std::max_element( variants.cbegin(), variants.cend(), endComp );
            const auto end = ( *maxElement )->end();
            const auto frontPadding = start - paddedRegion.start();
            const auto backPadding = paddedRegion.end() - end;

            region = caller::Region( paddedRegion.contig(), utils::Interval( start, end ) );
            referenceString = m_referenceSequence->subseq( region ).sequence();
            altString = sequence.substr( frontPadding, sequence.size() - frontPadding - backPadding );
        }

        const std::size_t matrixSize = ( region.size() + 1 ) * ( altString.size() + 1 );
        if ( matrixSize > m_maxMatrixSize )
        {
            ECHIDNA_LOG( DEBUG, "Skipping normalization because the matrix has size "
                                    << region.size() + 1 << "x" << altString.size() + 1 << " = " << matrixSize
                                    << ", which is greater than the threshold of " << m_maxMatrixSize );
            if ( unnormalizedVariants.is_initialized() )
            {
                return unnormalizedVariants.get();
            }
            else
            {
                return {};
            }
        }

        const auto nWPenalties = utils::NWPenalties();
        auto needlemanWunsch = utils::NeedlemanWunsch( referenceString, altString, nWPenalties );

        needlemanWunsch.getScoreMatrix();

        auto nWVariants = needlemanWunsch.traceBack();

        variantSet_t variantSet = {};

        for ( auto & simpleVariant : nWVariants )
        {
            variantSet.insert( simpleVariant.getVariant( region.contig(), region.start(), m_referenceSequence )
                               //                                      ->getLeftAligned( paddedRegion.start() )
                               );
        }

        return variantSet;
    }
}
}