// All content Copyright (C) 2018 Genomics plc
#ifndef GENOMICS_ECHIDNA_LOGGING_HPP
#define GENOMICS_ECHIDNA_LOGGING_HPP

#include "boost/log/trivial.hpp"
#include "boost/log/sources/logger.hpp"
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/severity_channel_logger.hpp>

#include "common.hpp"
#include "utils/exceptions.hpp"

namespace echidna
{
namespace utils
{
    static boost::log::sources::severity_logger_mt< loggingLevel > g_echidnaLog;

    /// Sets up the logging configuration for Echidna
    ///
    /// @param loggingParams Logging parameter list from Echidna
    void initialiseLog( int logLevel, std::string logFileName, bool quietMode, int verbosity, bool logTimings );
}
}

#define ECHIDNA_LOG( LEVEL, MESSAGE )                                         \
    do                                                                        \
    {                                                                         \
        BOOST_LOG_SEV( utils::g_echidnaLog, loggingLevel::LEVEL ) << MESSAGE; \
    } while ( 0 )

#define ECHIDNA_ASSERT( COND, MESSAGE )                \
    do                                                 \
    {                                                  \
        if ( not( COND ) )                             \
        {                                              \
            ECHIDNA_LOG( FATAL, MESSAGE );             \
            throw utils::echidna_exception( MESSAGE ); \
        }                                              \
    } while ( 0 )

#define ECHIDNA_ERROR( COND, MESSAGE )                 \
    do                                                 \
    {                                                  \
        if ( not( COND ) )                             \
        {                                              \
            ECHIDNA_LOG( ERROR, MESSAGE );             \
            throw utils::echidna_exception( MESSAGE ); \
        }                                              \
    } while ( 0 )

#endif
