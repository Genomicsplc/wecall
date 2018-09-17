// All content Copyright (C) 2018 Genomics plc
#ifndef VARIANT_GENERATION_DATA_HPP
#define VARIANT_GENERATION_DATA_HPP

#include <string>
#include <memory>
#include <cstdint>
#include <utils/sequence.hpp>
#include <utils/referenceSequence.hpp>

#include "caller/region.hpp"
#include "common.hpp"

namespace echidna
{
namespace io
{
    class Read;
}  // namespace  io

namespace variant
{
    struct VariantGenerationData
    {
        explicit VariantGenerationData( utils::referenceSequencePtr_t refSeq,
                                        int64_t readStartPos,
                                        utils::BasePairSequence readSeq )
            : refSeq( refSeq ), readStartPos( readStartPos ), readSeq( readSeq )
        {
        }

        VariantGenerationData( utils::referenceSequencePtr_t refSeq, const io::Read & read );

        utils::referenceSequencePtr_t refSeq;

        const int64_t readStartPos;
        const utils::BasePairSequence readSeq;
    };

    using variantGenerationDataPtr_t = std::shared_ptr< const VariantGenerationData >;

}  // namespace variant
}  // namespace echidna

#endif  // VARIANT_GENERATION_DATA_HPP
