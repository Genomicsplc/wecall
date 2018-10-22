// All content Copyright (C) 2018 Genomics plc
#include "variantGenerationData.hpp"
#include "io/read.hpp"

namespace wecall
{
namespace variant
{
    VariantGenerationData::VariantGenerationData( utils::referenceSequencePtr_t refSeq, const io::Read & read )
        : VariantGenerationData( refSeq, read.getStartPos(), read.sequence() )
    {
    }

    //-----------------------------------------------------------------------------------------

}  // namespace variant
}  // namespace wecall
