// All content Copyright (C) 2018 Genomics plc
#ifndef WECALL_VARIANTCOMBINATIONS_H
#define WECALL_VARIANTCOMBINATIONS_H

#include "variant/type/variant.hpp"

namespace echidna
{
namespace variant
{
    enum VariantPairCombinationState
    {
        ALWAYS_TOGETHER,
        NEVER_TOGETHER,
        FIRST_IMPLIES_SECOND,
        SECOND_IMPLIES_FIRST,
        UNCERTAIN
    };

    class VariantCombinations
    {
    public:
        VariantCombinations( std::size_t minReadsToSupportClaim, int64_t maxClusterDistance );

        std::size_t nVariantCombinations() const { return m_variantCombinations.size(); }
        const std::vector< std::vector< varPtr_t > > & variantCombinations() const { return m_variantCombinations; }
        void setVariantCombinations( const std::vector< std::vector< varPtr_t > > & combinations )
        {
            m_variantCombinations = combinations;
        }
        bool allCombinationsComputed() const { return m_allCombinationsComputed; }

        void computeVariantCombinations( const std::vector< varPtr_t > & variants,
                                         const size_t maxCombinations,
                                         const utils::referenceSequencePtr_t & reference );

        std::vector< std::vector< varPtr_t > > filterVariantCombinations(
            const std::vector< std::vector< varPtr_t > > & variantCombinations,
            const std::vector< varPtr_t > & alwaysTogetherVariants,
            const std::vector< varPtr_t > & neverTogetherVariants,
            varPtr_t currentVariant );

        VariantPairCombinationState getState( varPtr_t first, varPtr_t second ) const;

        void getOverlappingReads( const varPtr_t & first,
                                  const varPtr_t & second,
                                  std::vector< io::readPtr_t > & firstReads,
                                  std::vector< io::readPtr_t > & secondReads ) const;

        std::size_t minReadsToSupportClaim() const { return m_minReadsToSupportClaim; }
        int64_t maxClusterDistance() const { return m_maxClusterDistance; }

    private:
        std::vector< std::vector< varPtr_t > > m_variantCombinations;

        bool m_allCombinationsComputed = false;

        std::size_t m_minReadsToSupportClaim;
        std::size_t m_maxClusterDistance;
    };
}
}

#endif  // WECALL_VARIANTCOMBINATIONS_H
