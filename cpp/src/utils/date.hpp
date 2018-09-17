// All content Copyright (C) 2018 Genomics plc
#ifndef UTILS_DATE_HPP
#define UTILS_DATE_HPP

#include "common.hpp"

#include <ctime>

namespace echidna
{
namespace utils
{
    /// Utility function to return the current date as ISO8601 formated string
    /// TODO: replace this messy C implementation with something nicer, maybe BOOST?
    inline std::string getCurrentDate()
    {
        char buffer[12];

        auto time = std::time( nullptr );
        auto tm = std::localtime( &time );
        std::strftime( buffer, 12, "%Y-%m-%d", tm );

        std::string date = buffer;
        return date;
    }
}
}
#endif
