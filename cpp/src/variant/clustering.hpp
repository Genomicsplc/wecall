// All content Copyright (C) 2018 Genomics plc
#ifndef CLUSTERING_HPP
#define CLUSTERING_HPP

#include <cmath>
#include "common.hpp"
#include "variant/type/variant.hpp"
#include "caller/region.hpp"
#include "caller/callSet.hpp"
#include "variant/variantContainer.hpp"
#include "utils/logging.hpp"
#include "variantCombinations.hpp"

namespace echidna
{
namespace variant
{
    class VariantCluster
    {
    public:
        VariantCluster();
        VariantCluster( const std::vector< varPtr_t > & variants, const caller::SetRegions & region );
        VariantCluster & operator+=( const VariantCluster & rhs );

        void push_back( const varPtr_t & varPtr, const caller::SetRegions & variantRegions );
        bool empty() const { return m_variants.empty(); }

        const std::vector< varPtr_t > & variants() const { return m_variants; };
        const caller::SetRegions & readRegions() const { return m_variantRegions; }
        caller::Region region() const { return m_variantRegions.getSpan(); }

        void setPaddedRegion( const caller::Region & paddedRegion ) { m_paddedRegion = paddedRegion; }
        const caller::Region & paddedRegion() const { return m_paddedRegion; }

        std::string toString() const;
        void sort() { std::sort( m_variants.begin(), m_variants.end(), varPtrComp() ); }

        void expandRegions( const int64_t desiredPadding );

        std::tuple< VariantCluster, std::vector< VariantCluster > > buildSubClusters(
            const int64_t clusterDistance ) const;

        int64_t zeroIndexedVcfStart( caller::callVector_t calls ) const
        {
            int64_t clusterStart = this->region().start();
            if ( m_zeroIndexedVcfPosition != this->region().start() )
            {
                auto firstVarPos = caller::zeroIndexedVCFPosition( calls );
                if ( m_zeroIndexedVcfPosition == firstVarPos )
                {
                    clusterStart = m_zeroIndexedVcfPosition;
                }
            }
            return clusterStart;
        }

        void computeVariantCombinations( std::size_t minReadsToSupportClaim,
                                         int64_t correlatedDistance,
                                         const std::size_t maxCombinations,
                                         const utils::referenceSequencePtr_t & reference )
        {
            VariantCombinations combsGen( minReadsToSupportClaim, correlatedDistance );
            combsGen.computeVariantCombinations( this->variants(), maxCombinations, reference );

            this->setVariantCombinations( combsGen.variantCombinations(), combsGen.allCombinationsComputed() );
        }

        void setVariantCombinations( const std::vector< std::vector< varPtr_t > > & variantCombinations,
                                     bool allCombinationsComputed )
        {
            m_variantCombinations = variantCombinations;
            m_allCombinationsComputed = allCombinationsComputed;
        }
        std::size_t nVariantCombinations() const { return m_variantCombinations.size(); }
        const std::vector< std::vector< varPtr_t > > & variantCombinations() const { return m_variantCombinations; }
        bool allCombinationsComputed() const { return m_allCombinationsComputed; }

    private:
        std::vector< varPtr_t > m_variants;
        caller::SetRegions m_variantRegions;
        caller::Region m_paddedRegion;
        int64_t m_zeroIndexedVcfPosition;

        std::vector< std::vector< varPtr_t > > m_variantCombinations;
        bool m_allCombinationsComputed;
    };

    void removeClustersNearBlockEnd( caller::Region & blockRegion,
                                     std::vector< VariantCluster > & clusters,
                                     const int maxClusterDist );

    std::vector< caller::Region > computeClustersPaddingRegions( const caller::Region & blockRegion,
                                                                 const std::vector< VariantCluster > & clusters );

    std::vector< VariantCluster > mergeClusters( const std::vector< VariantCluster > & input,
                                                 const int64_t maxVariantCombinations,
                                                 const int64_t maxDistBetweenClusters,
                                                 const int64_t maxClusterSize );

    std::vector< VariantCluster > generateVariantClusters( const variantSet_t & variants,
                                                           const int64_t minDistance,
                                                           const caller::Region & region );

    std::vector< VariantCluster > generateMergedClusters( const variantSet_t & variants,
                                                          const int64_t minDistance,
                                                          const int64_t maxDistance,
                                                          const int64_t maxClusterSize,
                                                          const caller::Region & region,
                                                          utils::referenceSequencePtr_t referenceSequence,
                                                          const std::size_t maxCombinations,
                                                          const std::size_t minReadsToSupportClaim );
}
}

#endif
