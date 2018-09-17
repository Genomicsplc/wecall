// All content Copyright (C) 2018 Genomics plc
#ifndef VARIANT_GENERATOR_HPP
#define VARIANT_GENERATOR_HPP

#include "caller/region.hpp"
#include "io/readRange.hpp"
#include "variant/type/variant.hpp"
#include "variant/variantContainer.hpp"
#include "utils/referenceSequence.hpp"

#include <vector>
#include <string>

namespace echidna
{
namespace variant
{
    std::vector< varPtr_t > normaliseVariantsOnStrand( const std::vector< varPtr_t > & variants,
                                                       const utils::ReferenceSequence & m_referenceSequence );

    std::vector< phred_t > getVariantsReadBaseQualities( int64_t startPos,
                                                         const utils::QualitySequence & qualityString,
                                                         const std::vector< varPtr_t > & variants,
                                                         const std::vector< variant::breakpointPtr_t > & breakpoints );

    class VariantGenerator
    {
    public:
        VariantGenerator( const utils::referenceSequencePtr_t & refSeq,
                          const phred_t minBaseQual,
                          const phred_t minMappingQual );

        VariantContainer generateVariantsFromReads( io::perSampleRegionsReads_t perSamReadRanges ) const;

    private:
        const utils::referenceSequencePtr_t m_referenceSequence;
        const phred_t m_minBaseQual;
        const phred_t m_minMappingQual;
    };
}
}

#endif
