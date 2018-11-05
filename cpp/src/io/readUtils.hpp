// All content Copyright (C) 2018 Genomics plc
#ifndef READ_UTILS_HPP
#define READ_UTILS_HPP

#include "common.hpp"
#include "io/read.hpp"

namespace wecall
{
namespace io
{
    namespace read
    {
        phred_t minBaseQualityInReadAroundInterval( const Read & read,
                                                    const utils::Interval & interval,
                                                    const int64_t padding );
    }
}
}
#endif