// All content Copyright (C) 2018 Genomics plc
#include "utils/logging.hpp"

#include <fstream>

#include <boost/log/common.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/attributes.hpp>
#include <boost/log/sinks.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/support/date_time.hpp>

#include "caller/params.hpp"

namespace echidna
{
namespace utils
{
    namespace logging = boost::log;
    namespace expr = boost::log::expressions;
    namespace sinks = boost::log::sinks;
    namespace keywords = boost::log::keywords;
    namespace src = boost::log::sources;

    BOOST_LOG_ATTRIBUTE_KEYWORD( severity, "Severity", loggingLevel )

    struct empty_deleter
    {
        //! Function object result type
        typedef void result_type;
        /*!
         * Does nothing
         */
        void operator()( const volatile void * ) const {}
    };

    void initialiseLog( int logLevel, std::string logFileName, bool quietMode, int verbosity, bool logTimings )
    {
        if ( logTimings and logLevel < loggingLevel::TIMING )
        {
            logLevel = loggingLevel::TIMING;
        }
        if ( quietMode )
        {
            verbosity = loggingLevel::FATAL;
        }
        // TODO: Make sure that everything works when run with multiple threads.

        // Set up logging output for Echidna. The various output handlers are added
        // to the logging core,
        using sink_t = sinks::synchronous_sink< sinks::text_ostream_backend >;

        // The core singleton
        auto core = logging::core::get();

        // Create a simple clog-style stream backend
        auto backend = boost::make_shared< sinks::text_ostream_backend >();
        backend->add_stream( boost::shared_ptr< std::ostream >( &std::clog, empty_deleter() ) );

        // Flush stream after each log record is written.
        backend->auto_flush( true );
        auto sink( boost::make_shared< sink_t >( backend ) );

        sink->set_filter( severity <= verbosity );

        sink->set_formatter( expr::stream
                             << expr::format_date_time< boost::posix_time::ptime >( "TimeStamp", "%Y-%m-%d %H:%M:%S" )
                             << " -- " << severity << " -- " << expr::smessage );

        core->add_sink( sink );

        // Create a log file backend. TODO: check if std::ostream is thread-safe and replace
        // with the boost text_file_backend if necessary
        auto backend2 = boost::make_shared< sinks::text_ostream_backend >();
        backend2->add_stream( boost::shared_ptr< std::ostream >( new std::ofstream( logFileName ) ) );
        backend2->auto_flush( true );
        auto sink2( boost::make_shared< sink_t >( backend2 ) );
        sink2->set_filter( severity <= logLevel );

        sink2->set_formatter( expr::stream
                              << expr::format_date_time< boost::posix_time::ptime >( "TimeStamp", "%Y-%m-%d %H:%M:%S" )
                              << " -- " << severity << " -- " << expr::smessage );

        core->add_sink( sink2 );
        logging::add_common_attributes();
    }
}
}
