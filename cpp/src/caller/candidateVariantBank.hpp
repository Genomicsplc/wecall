// All content Copyright (C) 2018 Genomics plc
#ifndef CANDIDATE_VARIANT_BANK_HPP
#define CANDIDATE_VARIANT_BANK_HPP

#include <string>

#include "caller/region.hpp"
#include "utils/logging.hpp"
#include "variant/type/variant.hpp"
#include "vcf/record.hpp"
#include "io/fastaFile.hpp"

namespace wecall
{
namespace caller
{
    std::vector< double > getPriorsFromInfo( const vcf::Info & info );
}
}

#endif
