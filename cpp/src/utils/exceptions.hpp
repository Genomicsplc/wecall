// All content Copyright (C) 2018 Genomics plc
#ifndef UTILS_EXCEPTIONS_HPP
#define UTILS_EXCEPTIONS_HPP

#include <stdexcept>

namespace wecall
{
namespace utils
{
    class wecall_exception : public std::runtime_error
    {
    public:
        wecall_exception( const std::string & error ) : std::runtime_error( error ) {}
    };
}
}

#endif
