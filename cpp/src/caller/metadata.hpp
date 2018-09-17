// All content Copyright (C) 2018 Genomics plc
#ifndef METADATA_HPP
#define METADATA_HPP
#include "common.hpp"

namespace echidna
{
namespace caller
{
    struct GenotypeMetadata
    {
        phred_t phaseQuality = constants::unknownValue;
        phred_t genotypeQuality = constants::unknownValue;
    };
}
}

#endif
