// All content Copyright (C) 2018 Genomics plc
#ifndef UTILS_EXCEPTIONS_HPP
#define UTILS_EXCEPTIONS_HPP

#include <stdexcept>

namespace echidna
{
namespace utils
{
    class echidna_exception : public std::runtime_error
    {
    public:
        echidna_exception( const std::string & error ) : std::runtime_error( error ) {}
    };
}
}

#endif
