// All content Copyright (C) 2018 Genomics plc
#ifndef VARIANT_SNPFINDER_HPP
#define VARIANT_SNPFINDER_HPP

#include "variant/variantGenerationData.hpp"
#include "alignment/cigarItems.hpp"
#include "variant/type/variant.hpp"

namespace wecall
{
namespace variant
{
    class SNPFinder
    {
    public:
        SNPFinder( variantGenerationDataPtr_t variantGenerationData ) : m_varGenData( variantGenerationData ) {}

        variantSet_t findSNPsInReadSegment( const alignment::offsetsPtr_t offsets, const int64_t length ) const;

    private:
        int64_t refIndexFromReadIndex( const int64_t readIndex, const alignment::offsetsPtr_t offsetsPtr_t ) const;
        int64_t readIndexFromRefIndexInStr( const alignment::offsetsPtr_t offsetsPtr_t ) const;

        variantGenerationDataPtr_t m_varGenData;
    };
}
}

#endif
