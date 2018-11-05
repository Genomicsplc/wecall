// All content Copyright (C) 2018 Genomics plc
#include <string>
#include "utils/exceptions.hpp"

namespace wecall
{
namespace test
{
    std::string requireEnv( const char * variableName )
    {
        char * value = std::getenv( variableName );
        if ( value == nullptr )
        {
            throw utils::wecall_exception( std::string( "Undefined environment variable: " ) + variableName );
        }
        else
        {
            return std::string( value );
        }
    }
}
}
