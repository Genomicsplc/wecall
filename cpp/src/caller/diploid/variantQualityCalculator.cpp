// All content Copyright (C) 2018 Genomics plc
#include "caller/diploid/variantQualityCalculator.hpp"
#include "utils/indexedProduct.hpp"

#include <numeric>

namespace echidna
{
namespace caller
{
    namespace model
    {
        //-----------------------------------------------------------------------------------------

        std::vector< double > VariantQualityCalculator::getQualities() const
        {
            std::vector< double > variantQualities( m_candidateVariants.size() );
            for ( std::size_t varIndex = 0; varIndex < m_candidateVariants.size(); ++varIndex )
            {
                variantQualities[varIndex] = this->computeVariantQuality( varIndex );
            }
            return variantQualities;
        }

        //-----------------------------------------------------------------------------------------

        double VariantQualityCalculator::computeVariantQuality( const std::size_t variantIndex ) const
        {
            const auto variant = m_candidateVariants[variantIndex];
            const auto prior = variant->prior();
            const auto nSamples = m_samples.size();

            auto sumLogProbNoVariant = 0.0;
            auto sumLogProbTotal = 0.0;

            for ( std::size_t sampleIndex = 0; sampleIndex < nSamples; ++sampleIndex )
            {
                sumLogProbNoVariant += m_reweightedVariantQualities[sampleIndex][variantIndex];
                sumLogProbTotal += m_variantQualitiesTotal[sampleIndex];
            }

            const auto weightedRatio = std::max( std::numeric_limits< double >::min(),
                                                 exp( sumLogProbNoVariant - sumLogProbTotal ) * ( 1.0 - prior ) );
            const auto posterior =
                std::min( round( -10.0 * ( log10( weightedRatio ) - log10( prior + weightedRatio ) ) ),
                          constants::maxPhredScore );

            if ( false )
            {
                ECHIDNA_LOG( SUPER_DEBUG, "sumLogProbNoVariant " << sumLogProbNoVariant );
                ECHIDNA_LOG( SUPER_DEBUG, "sumLogTotalEvents " << sumLogProbTotal );
                ECHIDNA_LOG( SUPER_DEBUG, "Prior is " << prior << " for var " << variant->toString() );
                ECHIDNA_LOG( SUPER_DEBUG, "Ratio is " << weightedRatio << " for var " << variant->toString() );
                ECHIDNA_LOG( SUPER_DEBUG, "Posterior is " << posterior << " for var " << variant->toString() );
            }

            return posterior;
        }

        //-----------------------------------------------------------------------------------------
    }
}
}
