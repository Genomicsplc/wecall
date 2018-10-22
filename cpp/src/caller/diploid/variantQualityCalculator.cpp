// All content Copyright (C) 2018 Genomics plc
#include "caller/diploid/variantQualityCalculator.hpp"
#include "utils/indexedProduct.hpp"

#include <numeric>

namespace wecall
{
namespace caller
{
    namespace model
    {
        //-----------------------------------------------------------------------------------------

        std::vector< double > VariantQualityCalculator::computeReweightedHaplotypeFrequencies(
            const std::size_t variantIndex ) const
        {
            std::vector< double > reweightedFrequenciesWithoutVariant = m_haplotypeFrequencies;

            const auto hapIndexesThisVariant =
                m_candidateHaplotypes.getHaplotypeIndicesForVariant( m_candidateVariants[variantIndex] );

            for ( const auto hapIndex : hapIndexesThisVariant )
            {
                reweightedFrequenciesWithoutVariant[hapIndex] = 0.0;
            }

            const auto sumFreqsWithoutVariant = std::accumulate( reweightedFrequenciesWithoutVariant.begin(),
                                                                 reweightedFrequenciesWithoutVariant.end(), 0.0 );

            if ( sumFreqsWithoutVariant > 0.0 )
            {
                for ( std::size_t haplotypeIndex = 0; haplotypeIndex < m_candidateHaplotypes.size(); ++haplotypeIndex )
                {
                    reweightedFrequenciesWithoutVariant[haplotypeIndex] /= sumFreqsWithoutVariant;
                }
            }

            return reweightedFrequenciesWithoutVariant;
        }

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
            const auto logOfMinDouble = log( std::numeric_limits< double >::min() );
            const auto & samples = m_samples;
            const auto nSamples = m_samples.size();

            const auto & reads = m_readRangesPerSample;

            const auto reweightedFrequenciesWithoutVariant =
                this->computeReweightedHaplotypeFrequencies( variantIndex );

            auto sumLogProbNoVariant = 0.0;
            auto sumLogTotalEvents = 0.0;

            for ( std::size_t sampleIndex = 0; sampleIndex < nSamples; ++sampleIndex )
            {
                const auto & sample = samples[sampleIndex];
                const auto nGenotypes = m_candidateGenotypes.at( sample ).size();
                const auto & readItPair = reads.at( sample );
                const auto nReads = std::distance( readItPair.begin(), readItPair.end() );

                if ( nReads == 0 )
                {
                    continue;
                }

                auto sumProbNoVariantThisIndividual = 0.0;
                auto sumProbTotalThisIndividual = 0.0;

                for ( std::size_t genotypeIndex = 0; genotypeIndex < nGenotypes; ++genotypeIndex )
                {
                    const auto hapIndicies = m_candidateGenotypes.at( sample ).getHaplotypeIndices( genotypeIndex );

                    const auto adjustedGenotypeLikelihood =
                        m_genotypeLikelihoods[sampleIndex][genotypeIndex] *
                        m_candidateGenotypes.at( sample )[genotypeIndex]->nCombinationsThisGenotype();

                    sumProbTotalThisIndividual += utils::indexedProduct(
                        m_haplotypeFrequencies, hapIndicies.cbegin(), hapIndicies.cend(), adjustedGenotypeLikelihood );
                    sumProbNoVariantThisIndividual +=
                        utils::indexedProduct( reweightedFrequenciesWithoutVariant, hapIndicies.cbegin(),
                                               hapIndicies.cend(), adjustedGenotypeLikelihood );
                }

                sumLogTotalEvents +=
                    ( sumProbTotalThisIndividual > 0 ) ? log( sumProbTotalThisIndividual ) : logOfMinDouble;
                sumLogProbNoVariant +=
                    ( sumProbNoVariantThisIndividual > 0 ) ? log( sumProbNoVariantThisIndividual ) : logOfMinDouble;
            }

            const auto weightedRatio = std::max( std::numeric_limits< double >::min(),
                                                 exp( sumLogProbNoVariant - sumLogTotalEvents ) * ( 1.0 - prior ) );
            const auto posterior =
                std::min( round( -10.0 * ( log10( weightedRatio ) - log10( prior + weightedRatio ) ) ),
                          constants::maxPhredScore );

            if ( false )
            {
                WECALL_LOG( SUPER_DEBUG, "sumLogProbNoVariant " << sumLogProbNoVariant );
                WECALL_LOG( SUPER_DEBUG, "sumLogTotalEvents " << sumLogTotalEvents );
                WECALL_LOG( SUPER_DEBUG, "Prior is " << prior << " for var " << variant->toString() );
                WECALL_LOG( SUPER_DEBUG, "Ratio is " << weightedRatio << " for var " << variant->toString() );
                WECALL_LOG( SUPER_DEBUG, "Posterior is " << posterior << " for var " << variant->toString() );
            }

            return posterior;
        }

        //-----------------------------------------------------------------------------------------
    }
}
}
