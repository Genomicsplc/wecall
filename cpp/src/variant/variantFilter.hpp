// All content Copyright (C) 2018 Genomics plc
#ifndef VARIANT_FILTER_HPP
#define VARIANT_FILTER_HPP

#include "variant/type/variant.hpp"

namespace wecall
{
namespace variant
{
    class VariantContainer;

    /// VariantFilter class. Class that from a variant container outputs sorted filtered variants
    ///
    class VariantFilter
    {

    public:
        /// Construct a VariantFilter class
        VariantFilter( const int64_t minReads, const int64_t minPerSamplePercentage )
            : m_minReads( minReads ), m_minPerSamplePercentage( minPerSamplePercentage )
        {
        }

        /// Returns sorted vector of all variants that pass filtering
        ///
        /// @param  varContainer The input variants to be filtered and sorted.
        variantSet_t getSortedFilteredVariants( const caller::Region & blockRegion,
                                                const VariantContainer & varContainer ) const;

    private:
        /// Returns true if variant has >= minReads supporting it
        bool variantPassesFilters( varPtr_t varPtr, const VariantContainer & varContainer ) const;

        bool variantHasMinReadsAcrossSamples( varPtr_t varPtr, const VariantContainer & varContainer ) const;

        bool variantHasSampleWithMinPercentReadCoverage( varPtr_t varPtr, const VariantContainer & varContainer ) const;

    private:
        const int64_t m_minReads;
        const int64_t m_minPerSamplePercentage;
    };
}
}

#endif
