// All content Copyright (C) 2018 Genomics plc
#include <utils/bestScoreSelector.hpp>
#include <caller/params.hpp>
#include "caller/haplotypeLikelihoods.hpp"
#include "utils/partition.hpp"
#include "variant/haplotypeGenerator.hpp"
#include "variant/haplotype.hpp"
#include "utils/combinations.hpp"
#include "variant/haplotypeRanker.hpp"
#include "variant/type/variant.hpp"

namespace echidna
{
namespace variant
{

    AlignmentHaplotypeGenerator::AlignmentHaplotypeGenerator( const std::vector< variant::varPtr_t > & variants,
                                                              const caller::SetRegions & regions,
                                                              const io::perSampleRegionsReads_t & readsPerSample,
                                                              utils::referenceSequencePtr_t referenceSequence,
                                                              const int64_t maxHaplotypesPerRanker,
                                                              const std::size_t minReadsToSupportClaim )
        : m_vars( variants ),
          m_regions( regions ),
          m_readsPerSample( readsPerSample ),
          m_referenceSequence( referenceSequence ),
          m_maxHaplotypesPerRanker( maxHaplotypesPerRanker ),
          m_minReadsToSupportClaim( minReadsToSupportClaim )
    {
    }

    std::vector< variantSet_t > AlignmentHaplotypeGenerator::bestVariantCombos( const VariantCluster & cluster ) const
    {
        ECHIDNA_LOG( DEBUG, "AlignmentHaplotypeGenerator::bestVariantCombos for " << cluster.toString() );
        ECHIDNA_ASSERT( cluster.allCombinationsComputed(), "This step requires all combinations to be computed" );

        const auto clusterReads = io::reduceRegionSet( m_readsPerSample, cluster.readRegions() );
        HaplotypeVector haplotypes( cluster.readRegions(), m_referenceSequence );
        for ( const auto & varCombo : cluster.variantCombinations() )
        {
            haplotypes.push_back( variant::variantSet_t( varCombo.cbegin(), varCombo.cend() ) );
        }
        haplotypes.sort();
        haplotypes.merge();

        const AlignmentHaplotypeRanker alignmentHaplotypeRanker( clusterReads );

        const auto bestIndicies = alignmentHaplotypeRanker.getTopHaplotypes( haplotypes, m_maxHaplotypesPerRanker );

        std::vector< variantSet_t > goodVariantCombos;
        for ( const auto & goodIndex : bestIndicies )
        {
            goodVariantCombos.push_back( haplotypes[goodIndex].getVariants() );
        }
        return goodVariantCombos;
    }

    HaplotypeVector AlignmentHaplotypeGenerator::generateRawHaplotypes() const
    {
        HaplotypeVector haplotypes( m_regions, m_referenceSequence );

        ECHIDNA_LOG( DEBUG, "AlignmentHaplotypeGenerator::generateRawHaplotypes for " << m_regions );
        const auto region = m_regions.getSpan();

        variantSet_t setVariants( m_vars.cbegin(), m_vars.cend() );
        const auto reClusteredVariants = variant::generateMergedClusters(
            setVariants, constants::needlemanWunschPadding, region.size(), region.size(), region, m_referenceSequence,
            caller::params::defaults::maxClusterVariantCombinations, m_minReadsToSupportClaim );

        const auto exitPred = []( const VariantCluster & cluster )
        {
            return not cluster.allCombinationsComputed();
        };
        if ( std::any_of( reClusteredVariants.begin(), reClusteredVariants.end(), exitPred ) )
        {
            return haplotypes;
        }

        std::vector< std::vector< variantSet_t > > allVariantCombos;
        for ( const auto & miniCluster : reClusteredVariants )
        {
            const auto bestVariantCombos = this->bestVariantCombos( miniCluster );
            allVariantCombos.emplace_back( bestVariantCombos );
        }

        for ( const auto & newCombos : allVariantCombos )
        {
            std::vector< variantSet_t > currentCombos;
            for ( std::size_t hapIndex = 0; hapIndex < haplotypes.size(); ++hapIndex )
            {
                currentCombos.push_back( haplotypes[hapIndex].getVariants() );
            }
            for ( const auto & newCombo : newCombos )
            {
                for ( const auto & oldCombo : currentCombos )
                {
                    variantSet_t joinedCombo = oldCombo;
                    joinedCombo.insert( newCombo.cbegin(), newCombo.cend() );
                    if ( isValidVariantCombination( joinedCombo.cbegin(), joinedCombo.cend(), m_referenceSequence ) )
                    {
                        haplotypes.push_back( joinedCombo );
                    }
                }
                haplotypes.push_back( newCombo );
            }
            if ( haplotypes.size() >= static_cast< std::size_t >( m_maxHaplotypesPerRanker - 1 ) )
            {
                const AlignmentHaplotypeRanker alignmentHaplotypeRanker( m_readsPerSample );
                const auto bestIndicies =
                    alignmentHaplotypeRanker.getTopHaplotypes( haplotypes, m_maxHaplotypesPerRanker - 1 );
                haplotypes.keepIndicies( bestIndicies );
            }
        }

        return haplotypes;
    }

    HaplotypeVector AlignmentHaplotypeGenerator::generateHaplotypes() const
    {
        HaplotypeVector haplotypes = this->generateRawHaplotypes();

        const variantSet_t referenceHaplotypeVariantSet = {};
        haplotypes.push_back( referenceHaplotypeVariantSet );

        haplotypes.sort();
        haplotypes.merge();

        ECHIDNA_LOG( DEBUG, "Generated " << haplotypes.size() << " haplotypes." );
        ECHIDNA_LOG( SUPER_DEBUG, "Generated haplotype candidates:-" << std::endl
                                                                     << haplotypes.toString() );
        return haplotypes;
    }

    //-------------------------------------------------------------------------------------------------
}
}
