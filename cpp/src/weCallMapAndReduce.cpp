// All content Copyright (C) 2018 Genomics plc
#include "weCallMapAndReduce.hpp"

namespace echidna
{
using namespace boost::program_options;
using namespace io;

weCallMapAndReduce::weCallMapAndReduce()
    : weCallBase( constants::weCallString ),
      m_configOpts( "Config parameters" ),
      m_cmdLineOpts( "Command line parameters" )
{
    initOptions();
}

int weCallMapAndReduce::processJob( int argc, char * argv[] )
{
    try
    {
        // Get processing options
        variables_map optValues;
        caller::params::VcfOptionsLog vcfOptionsLog;

        // Note: Command line options, being first, override config file options
        basic_parsed_options< char > cmdLineValues =
            command_line_parser( argc, argv )
                .options( m_cmdLineOpts )
                .style( command_line_style::default_style & ~command_line_style::allow_guessing )
                .run();
        vcfOptionsLog.addOptions( cmdLineValues );
        store( cmdLineValues, optValues );

        if ( optValues.count( "help" ) )
        {
            std::cout << m_publicOpts << std::endl;
            return 0;
        }

        if ( optValues.count( "config" ) )
        {
            std::string configFileName = optValues["config"].as< std::string >();

            std::ifstream configIn( configFileName );
            if ( not configIn.is_open() )
            {
                throw utils::echidna_exception( "Could not open configuration file: " + configFileName );
            }

            // Options from the config file are added to those from the command line (no override)
            basic_parsed_options< char > configFileValues = parse_config_file( configIn, m_configOpts );
            vcfOptionsLog.addOptions( configFileValues );
            store( configFileValues, optValues );
        }

        // Parsed values are checked for completeness here - exceptions may be thrown
        notify( optValues );

        // Convert options to job parameters
        caller::params::Logging loggingParams( optValues );
        utils::initialiseLog( loggingParams.m_logLevel, loggingParams.m_logFilename, loggingParams.m_quietMode,
                              loggingParams.m_verbosity, loggingParams.m_logTimings );

        caller::params::Application applicationParams( m_programName, m_programVersion, g_GIT_SHA1, g_BUILD_DATE,
                                                       vcfOptionsLog.str() );

        caller::params::PrivateSystem privateSystemParams( optValues );
        caller::params::Data dataParams( optValues, privateSystemParams.m_overwrite );
        caller::params::System systemParams( optValues );
        caller::params::Filters filterParams( optValues );
        caller::params::PrivateCalling privateCallingParams( optValues );
        caller::params::Calling callingParams( optValues );
        caller::params::PrivateData privateOutputParams( optValues );

        if ( systemParams.m_numberOfJobs == 0 )
        {
            ECHIDNA_LOG( INFO, "Will run in serial mode" );

            caller::Job job( applicationParams, dataParams, systemParams, privateSystemParams, filterParams,
                             privateCallingParams, callingParams, privateOutputParams );

            job.process();
        }
        else
        {
            const caller::params::Reduce reduceParams( dataParams.workDir(), dataParams.outputDataSink() );
            caller::params::validateReduceParamsPreMap( reduceParams );

            std::shared_ptr< boost::asio::io_service > ioService( new boost::asio::io_service );
            std::shared_ptr< boost::asio::io_service::work > work( new boost::asio::io_service::work( *ioService ) );

            ECHIDNA_LOG( INFO, "Will run " << systemParams.m_numberOfJobs << " jobs simultaneously, where possible" );
            boost::thread_group threadPool;
            for ( std::size_t i = 0; i < systemParams.m_numberOfJobs; ++i )
            {
                threadPool.create_thread( boost::bind( &boost::asio::io_service::run, ioService ) );
            }

            const std::vector< caller::params::Data > chunkedDataParams = dataParams.splitWorkload();

            for ( std::size_t i = 0; i < chunkedDataParams.size(); ++i )
            {
                const auto jobProcessor = [&]( std::size_t jobIndex ) -> void
                {
                    ECHIDNA_LOG( INFO, "Started job " << jobIndex );

                    caller::Job job( applicationParams, chunkedDataParams[jobIndex], systemParams, privateSystemParams,
                                     filterParams, privateCallingParams, callingParams, privateOutputParams );

                    job.process();

                    ECHIDNA_LOG( INFO, "Finished job " << jobIndex );
                };

                ioService->post( boost::lambda::bind( jobProcessor, i ) );

                ECHIDNA_LOG( INFO, "Posted job " << i );
            }

            work.reset();

            ECHIDNA_LOG( INFO, "Waiting for all jobs to finish" );

            threadPool.join_all();

            ECHIDNA_LOG( INFO, "Joined all threads" );

            caller::JobReduce( reduceParams ).process();
        }

        ECHIDNA_LOG( DEBUG, m_programName + " completed" );
    }
    catch ( std::exception & e )
    {
        std::cerr << "FAILED - " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

void weCallMapAndReduce::initOptions()
{
    // Note: Because there is no standard way of conveying the absence of an option to the working code, please
    // ensure
    //       that every option is either required, has a default_value or is multi-value (an empty vector conveys
    //       absence).
    //
    // Note: default_value is the value used regardless of whether or not the option token is present.
    //       implicit_value is the value used if the option token is present but no value is specified.
    //       The latter should be included for all bool options.
    //
    // Note: Multi-value parameters are read into strings and parsed internally due to the limitations of the
    //       boost multitoken facility, which does not support single line lists in config files or comma
    //       delimitation.

    options_description basicOpts( "Basic Parameters", caller::params::PARAM_HELP_DISPLAY_WIDTH );
    basicOpts.add_options()( "help,h", "produce help message" )( "config,c", value< std::string >(),
                                                                 "configuration file name" );

    m_publicOpts.add( basicOpts )
        .add( caller::params::Data::getOptionsDescription() )
        .add( caller::params::Calling::getOptionsDescription() )
        .add( caller::params::Logging::getOptionsDescription() )
        .add( caller::params::System::getOptionsDescription() );

    m_configOpts.add( caller::params::Data::getOptionsDescription() )
        .add( caller::params::Calling::getOptionsDescription() )
        .add( caller::params::Logging::getOptionsDescription() )
        .add( caller::params::System::getOptionsDescription() )
        .add( caller::params::PrivateSystem::getOptionsDescription() )
        .add( caller::params::PrivateData::getOptionsDescription() )
        .add( caller::params::Filters::getOptionsDescription() )
        .add( caller::params::PrivateCalling::getOptionsDescription() );

    m_cmdLineOpts.add( basicOpts ).add( m_configOpts );
}
}
