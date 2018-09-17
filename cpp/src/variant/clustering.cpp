// All content Copyright (C) 2018 Genomics plc
#include "variant/clustering.hpp"
#include "haplotype.hpp"
#include "variantCombinations.hpp"

namespace echidna
{
namespace variant
{
    VariantCluster::VariantCluster()
        : m_variants(),
          m_variantRegions(),
          m_paddedRegion( "", -1, -1 ),
          m_zeroIndexedVcfPosition( -1 ),
          m_variantCombinations(),
          m_allCombinationsComputed( false )
    {
    }

    void VariantCluster::push_back( const varPtr_t & varPtr, const caller::SetRegions & variantRegions )
    {
        ECHIDNA_ASSERT( variantRegions.size() > 0, "Variant regions must be non-empty" );
        const auto span = variantRegions.getSpan();

        ECHIDNA_ASSERT( span.start() <= varPtr->start(),
                        "Variant start region should be a padding of the chosen representation start pos." );
        ECHIDNA_ASSERT( span.end() >= varPtr->end(),
                        "Variant end region should be a padding of the chosen representation end pos." );
        const auto varZeroStart = std::min( span.start(), varPtr->zeroIndexedVcfPosition() );
        if ( m_variants.empty() )
        {
            m_variants = {varPtr};
            for ( const auto & reg : variantRegions )
            {
                m_variantRegions.insert( reg );
            }
            m_paddedRegion = m_variantRegions.getSpan();
            m_zeroIndexedVcfPosition = varZeroStart;
        }
        else
        {
            m_variants.push_back( varPtr );
            for ( const auto & reg : variantRegions )
            {
                m_variantRegions.insert( reg );
                m_paddedRegion.combine( reg );
            }

            m_zeroIndexedVcfPosition = std::min( m_zeroIndexedVcfPosition, varZeroStart );
        }
    }

    VariantCluster::VariantCluster( const std::vector< varPtr_t > & variants, const caller::SetRegions & regions )
        : m_variants(),
          m_variantRegions( regions ),
          m_paddedRegion( regions.getSpan() ),
          m_zeroIndexedVcfPosition( regions.getSpan().start() ),
          m_variantCombinations(),
          m_allCombinationsComputed( false )
    {

        for ( const auto & var : variants )
        {
            this->push_back( var, var->region() );
        }
    }

    VariantCluster & VariantCluster::operator+=( const VariantCluster & rhs )
    {
        ECHIDNA_ASSERT( this->region().end() <= rhs.region().start(), "Clusters should be in correct order." );

        this->m_variants.insert( this->m_variants.end(), rhs.m_variants.begin(), rhs.m_variants.end() );
        for ( const auto & region : rhs.m_variantRegions )
        {
            this->m_variantRegions.insert( region );
        }
        this->m_paddedRegion.combine( rhs.m_paddedRegion );
        this->setVariantCombinations( {}, false );
        return *this;
    }

    std::tuple< VariantCluster, std::vector< VariantCluster > > VariantCluster::buildSubClusters(
        const int64_t clusterDistance ) const
    {
        const auto thisRegion = region();
        const auto largeVariantThreshold = 2 * clusterDistance;
        const auto isLargeVar = [largeVariantThreshold]( const varPtr_t & var )
        {
            return var->sequenceLengthInRef() >= largeVariantThreshold;
        };
        const bool hasLargeVariant = std::any_of( this->variants().cbegin(), this->variants().cend(), isLargeVar );
        if ( not hasLargeVariant )
        {
            VariantCluster mainCluster = *this;
            std::vector< VariantCluster > clusters = {};
            return std::make_tuple( mainCluster, clusters );
        }

        variantSet_t smallVars;
        variantSet_t remaining;
        std::remove_copy_if( this->variants().cbegin(), this->variants().cend(),
                             std::inserter( smallVars, smallVars.end() ), isLargeVar );
        std::copy_if( this->variants().cbegin(), this->variants().cend(), std::inserter( remaining, remaining.end() ),
                      isLargeVar );

        auto smallClusters = generateVariantClusters( smallVars, clusterDistance, thisRegion );
        const auto paddingRegions = computeClustersPaddingRegions( thisRegion, smallClusters );
        for ( std::size_t clusterIndex = 0; clusterIndex < smallClusters.size(); ++clusterIndex )
        {
            smallClusters[clusterIndex].setPaddedRegion( paddingRegions[clusterIndex] );
        }

        variant::VariantCluster lvcCluster;
        for ( const auto & largeVariant : remaining )
        {
            const auto varRegions = largeVariant->getStartEndRegions( thisRegion.start(), thisRegion.end() );
            lvcCluster.push_back( largeVariant, varRegions );
        }

        std::vector< variant::VariantCluster > smallVariantClustersNotTouchingLargeVariants;
        for ( const auto & smallCluster : smallClusters )
        {
            if ( lvcCluster.readRegions().overlaps( smallCluster.region().getPadded( clusterDistance + 1 ) ) )
            {
                for ( const auto & smallVar : smallCluster.variants() )
                {
                    lvcCluster.push_back( smallVar, smallCluster.region() );
                }
            }
            else
            {
                smallVariantClustersNotTouchingLargeVariants.push_back( smallCluster );
            }
        }

        lvcCluster.setPaddedRegion( this->paddedRegion() );
        lvcCluster.sort();

        return std::make_tuple( lvcCluster, smallVariantClustersNotTouchingLargeVariants );
    }

    void VariantCluster::expandRegions( const int64_t desiredPadding )
    {
        m_variantRegions.fill( desiredPadding );

        const auto currentRegion = this->region();
        if ( m_paddedRegion.end() - desiredPadding > currentRegion.end() )
        {
            const auto newEnd = std::min( currentRegion.end() + desiredPadding, m_paddedRegion.end() - desiredPadding );
            m_variantRegions.insert( caller::Region( currentRegion.contig(), currentRegion.end(), newEnd ) );
        }
        if ( m_paddedRegion.start() + desiredPadding < currentRegion.start() )
        {
            const auto newStart =
                std::max( currentRegion.start() - desiredPadding, m_paddedRegion.start() + desiredPadding );
            m_variantRegions.insert( caller::Region( currentRegion.contig(), newStart, currentRegion.start() ) );
            m_zeroIndexedVcfPosition = std::min( m_zeroIndexedVcfPosition, region().start() );
        }
    }

    //-------------------------------------------------------------------------------------------------

    std::vector< VariantCluster > generateVariantClusters( const variantSet_t & variants,
                                                           const int64_t minDistance,
                                                           const caller::Region & region )
    {
        ECHIDNA_LOG( DEBUG, "Generating variant clusters for: " + region.toString() );
        /// Create and return a vector of VariantCluster s. Each VariantCluster is made up
        /// of variants which are either all at the same position, or all overlap, or are within
        /// a specified distance of each other.

        std::vector< VariantCluster > clusters;
        VariantCluster current;

        for ( const auto & var : variants )
        {
            const auto varRegions = var->getStartEndRegions( region.start(), region.end() );

            if ( current.empty() )
            {
                current.push_back( var, varRegions );
            }

            else
            {
                while ( not clusters.empty() and
                        ( varRegions.getSpan().start() - clusters.back().region().end() ) <= minDistance )
                {
                    clusters.back() += current;
                    current = clusters.back();
                    clusters.pop_back();
                }

                const auto distance = varRegions.getSpan().start() - current.region().end();
                if ( distance > minDistance )
                {
                    clusters.push_back( current );
                    current = VariantCluster();
                }

                current.push_back( var, varRegions );
            }
        }

        if ( not current.empty() )
        {
            clusters.push_back( current );
        }

        ECHIDNA_LOG( DEBUG, "Generated  " << clusters.size() << " clusters." );

        const auto paddedClusterRegs = variant::computeClustersPaddingRegions( region, clusters );
        ECHIDNA_ASSERT( paddedClusterRegs.size() == clusters.size(), "" );
        for ( std::size_t clusterIndex = 0; clusterIndex < clusters.size(); ++clusterIndex )
        {
            clusters[clusterIndex].setPaddedRegion( paddedClusterRegs[clusterIndex] );
        }
        return clusters;
    }

    void removeClustersNearBlockEnd( caller::Region & blockRegion,
                                     std::vector< VariantCluster > & clusters,
                                     const int maxClusterDist )
    {
        int64_t actualEndBlock = blockRegion.end();
        while ( clusters.size() > 1 and blockRegion.end() - clusters.back().region().end() <= maxClusterDist )
        {
            actualEndBlock = clusters.back().region().start();
            ECHIDNA_LOG( DEBUG, "Cluster " << clusters.back().toString() << " will be considered in next block." );
            clusters.pop_back();
        }

        blockRegion = caller::Region( blockRegion.contig(), blockRegion.start(), actualEndBlock );
    }

    std::vector< caller::Region > computeClustersPaddingRegions( const caller::Region & blockRegion,
                                                                 const std::vector< VariantCluster > & clusters )
    {
        std::vector< caller::Region > paddedRegions;
        for ( std::size_t clusterIndex = 0u; clusterIndex < clusters.size(); ++clusterIndex )
        {
            paddedRegions.emplace_back(
                blockRegion.contig(),
                clusterIndex > 0u ? clusters[clusterIndex - 1].region().end() : blockRegion.start(),
                clusterIndex < clusters.size() - 1u ? clusters[clusterIndex + 1u].region().start()
                                                    : blockRegion.end() );
        }
        return paddedRegions;
    }
    //-----------------------------------------------------------------------------------------

    std::vector< VariantCluster > mergeClusters( const std::vector< VariantCluster > & input,
                                                 const int64_t maxVariantCombinations,
                                                 const int64_t maxDistBetweenClusters,
                                                 const int64_t maxClusterSize )
    {
        if ( input.empty() )
        {
            return {};
        }

        std::vector< VariantCluster > mergedClusters;
        VariantCluster currentCluster = input.front();

        for ( std::size_t clusterIndex = 1; clusterIndex < input.size(); ++clusterIndex )
        {
            const VariantCluster nextCluster( input[clusterIndex] );
            const int64_t nVariantCombinations =
                currentCluster.nVariantCombinations() * nextCluster.nVariantCombinations();
            const int64_t distance = nextCluster.region().start() - currentCluster.region().end();
            const int64_t predictedClusterSize = nextCluster.region().end() - currentCluster.region().start();

            if ( not( currentCluster.allCombinationsComputed() and nextCluster.allCombinationsComputed() ) )
            {
                mergedClusters.push_back( currentCluster );
                currentCluster = nextCluster;
            }
            else if ( nVariantCombinations > maxVariantCombinations or distance > maxDistBetweenClusters or
                      predictedClusterSize > maxClusterSize )
            {
                mergedClusters.push_back( currentCluster );
                currentCluster = nextCluster;
            }
            else
            {
                currentCluster += nextCluster;
            }
        }

        mergedClusters.push_back( currentCluster );

        return mergedClusters;
    }

    std::vector< VariantCluster > generateMergedClusters( const variantSet_t & variants,
                                                          const int64_t minDistance,
                                                          const int64_t maxDistance,
                                                          const int64_t maxClusterSize,
                                                          const caller::Region & region,
                                                          utils::referenceSequencePtr_t referenceSequence,
                                                          const std::size_t maxCombinations,
                                                          const std::size_t minReadsToSupportClaim )
    {
        auto clusters = variant::generateVariantClusters( variants, minDistance, region );

        for ( auto & cluster : clusters )
        {
            cluster.computeVariantCombinations( minReadsToSupportClaim, maxDistance, maxCombinations,
                                                referenceSequence );
        }
        auto mergeClusterDist = 2 * minDistance;
        while ( mergeClusterDist < maxDistance )
        {
            clusters = variant::mergeClusters( clusters, maxCombinations, mergeClusterDist, maxClusterSize );
            mergeClusterDist += minDistance;
            for ( auto & cluster : clusters )
            {
                cluster.computeVariantCombinations( minReadsToSupportClaim, maxDistance, maxCombinations,
                                                    referenceSequence );
            }
        }

        clusters = variant::mergeClusters( clusters, maxCombinations, maxDistance, maxClusterSize );

        for ( auto & cluster : clusters )
        {
            cluster.computeVariantCombinations( minReadsToSupportClaim, maxDistance, maxCombinations,
                                                referenceSequence );
        }
        for ( auto & cluster : clusters )
        {
            cluster.expandRegions( minDistance );
        }
        return clusters;
    }

    //-----------------------------------------------------------------------------------------

    std::string VariantCluster::toString() const
    {
        std::stringstream s;
        s << "VariantCluster(" << readRegions();
        s << " size=" << this->m_variants.size() << " {";

        for ( const auto & variant : this->m_variants )
        {
            s << variant->toString() << ", ";
        }
        s << "})";

        return s.str();
    }

    //-------------------------------------------------------------------------------------------------
}
}
