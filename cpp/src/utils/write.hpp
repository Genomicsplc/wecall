// All content Copyright (C) 2018 Genomics plc
#ifndef UTILS_WRITE_HPP
#define UTILS_WRITE_HPP

#include <sstream>
#include <iomanip>

#include "common.hpp"

namespace echidna
{
namespace utils
{

    template < typename Type >
    std::string toString( const Type & value )
    {
        std::stringstream ret;
        ret << std::setprecision( 9 ) << value;
        return ret.str();
    }
}
}
#endif
