// All content Copyright (C) 2018 Genomics plc
#include "weCallReduce.hpp"

namespace echidna
{
using namespace boost::program_options;
using namespace io;

weCallReduce::weCallReduce() : weCallBase( constants::weCallString + " reduce" ) { initOptions(); }

int weCallReduce::processJob( int argc, char * argv[] )
{
    try
    {
        // Get processing options
        variables_map optValues;

        // Note: Command line options, being first, override config file options
        basic_parsed_options< char > cmdLineValues =
            command_line_parser( argc, argv )
                .options( m_publicOpts )
                .style( command_line_style::default_style & ~command_line_style::allow_guessing )
                .run();

        store( cmdLineValues, optValues );

        if ( optValues.count( "help" ) )
        {
            std::cout << m_publicOpts << std::endl;
            return 0;
        }

        notify( optValues );
        // Convert options to job parameters
        caller::params::Logging loggingParams( optValues );
        utils::initialiseLog( loggingParams.m_logLevel, loggingParams.m_logFilename, loggingParams.m_quietMode,
                              loggingParams.m_verbosity, loggingParams.m_logTimings );

        caller::params::Reduce reduceParams( optValues );
        caller::params::validateReduceParams( reduceParams );
        caller::JobReduce( reduceParams ).process();
    }
    catch ( std::exception & e )
    {
        std::cerr << "FAILED - " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

void weCallReduce::initOptions()
{
    options_description basicOpts( "Basic Parameters", caller::params::PARAM_HELP_DISPLAY_WIDTH );
    basicOpts.add_options()( "help,h", "produce help message" );

    m_publicOpts.add( caller::params::Reduce::getOptionsDescription() );
    m_publicOpts.add( caller::params::Logging::getOptionsDescription() );
    m_publicOpts.add( basicOpts );
}
}
