// All content Copyright (C) 2018 Genomics plc
#ifndef GENOMICS_WECALL_LOGGING_HPP
#define GENOMICS_WECALL_LOGGING_HPP

#include "boost/log/trivial.hpp"
#include "boost/log/sources/logger.hpp"
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/severity_channel_logger.hpp>

#include "common.hpp"
#include "utils/exceptions.hpp"

namespace wecall
{
namespace utils
{
    static boost::log::sources::severity_logger_mt< loggingLevel > g_wecallLog;

    /// Sets up the logging configuration for weCall
    ///
    /// @param loggingParams Logging parameter list from weCall
    void initialiseLog( int logLevel, std::string logFileName, bool quietMode, int verbosity, bool logTimings );
}
}

#define WECALL_LOG( LEVEL, MESSAGE )                                         \
    do                                                                        \
    {                                                                         \
        BOOST_LOG_SEV( utils::g_wecallLog, loggingLevel::LEVEL ) << MESSAGE; \
    } while ( 0 )

#define WECALL_ASSERT( COND, MESSAGE )                \
    do                                                 \
    {                                                  \
        if ( not( COND ) )                             \
        {                                              \
            WECALL_LOG( FATAL, MESSAGE );             \
            throw utils::wecall_exception( MESSAGE ); \
        }                                              \
    } while ( 0 )

#define WECALL_ERROR( COND, MESSAGE )                 \
    do                                                 \
    {                                                  \
        if ( not( COND ) )                             \
        {                                              \
            WECALL_LOG( ERROR, MESSAGE );             \
            throw utils::wecall_exception( MESSAGE ); \
        }                                              \
    } while ( 0 )

#endif
