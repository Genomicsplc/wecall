// All content Copyright (C) 2018 Genomics plc
// test_main.cpp
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SmallVariantCaller
#include <cstdlib>
#include <string>
#include <boost/test/unit_test.hpp>
#include "caller/params.hpp"
#include "utils/logging.hpp"

struct LoggingGlobalFixture
{
    LoggingGlobalFixture()
    {
        std::string logFile = "unittest-cpp.log";
        std::cout << "Initialising log to file " << logFile << std::endl;
        echidna::caller::params::Logging loggingParams( loggingLevel::DEBUG, logFile, false, -1, false );
        echidna::utils::initialiseLog( loggingParams.m_logLevel, loggingParams.m_logFilename, loggingParams.m_quietMode,
                                       loggingParams.m_verbosity, loggingParams.m_logTimings );
    }
};

BOOST_GLOBAL_FIXTURE( LoggingGlobalFixture );
