// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>

#include "caller/params.hpp"

using namespace boost::program_options;

BOOST_AUTO_TEST_CASE( testGetParamListEmptyCase )
{
    variables_map variablesMap;
    notify( variablesMap );

    auto paramList = echidna::caller::params::getParamList( "TEST", variablesMap );
    BOOST_CHECK( paramList.empty() );
}

BOOST_AUTO_TEST_CASE( testGetParamListWithEmptyString )
{
    variables_map variablesMap;
    variable_value value( std::string( "" ), false );

    variablesMap.insert( std::make_pair( "TEST", value ) );
    notify( variablesMap );

    auto paramList = echidna::caller::params::getParamList( "TEST", variablesMap );
    std::vector< std::string > expectedParams = {};

    BOOST_CHECK_EQUAL_COLLECTIONS( expectedParams.begin(), expectedParams.end(), paramList.begin(), paramList.end() );
}

BOOST_AUTO_TEST_CASE( testGetParamListEmptyWithCommas )
{
    variables_map variablesMap;
    variable_value value( std::string( ",," ), false );

    variablesMap.insert( std::make_pair( "TEST", value ) );
    notify( variablesMap );

    auto paramList = echidna::caller::params::getParamList( "TEST", variablesMap );
    std::vector< std::string > expectedParams = {};

    BOOST_CHECK_EQUAL_COLLECTIONS( expectedParams.begin(), expectedParams.end(), paramList.begin(), paramList.end() );
}

BOOST_AUTO_TEST_CASE( testGetParamListCommaAtStartEnd )
{
    variables_map variablesMap;
    variable_value value( std::string( ",Hello," ), false );

    variablesMap.insert( std::make_pair( "TEST", value ) );
    notify( variablesMap );

    auto paramList = echidna::caller::params::getParamList( "TEST", variablesMap );
    std::vector< std::string > expectedParams = {"Hello"};

    BOOST_CHECK_EQUAL_COLLECTIONS( expectedParams.begin(), expectedParams.end(), paramList.begin(), paramList.end() );
}

BOOST_AUTO_TEST_CASE( testGetParamList )
{
    variables_map variablesMap;
    variable_value value( std::string( "He,, l,lo" ), false );

    variablesMap.insert( std::make_pair( "TEST", value ) );
    notify( variablesMap );

    auto paramList = echidna::caller::params::getParamList( "TEST", variablesMap );
    std::vector< std::string > expectedParams = {"He", "l", "lo"};

    BOOST_CHECK_EQUAL_COLLECTIONS( expectedParams.begin(), expectedParams.end(), paramList.begin(), paramList.end() );
}

BOOST_AUTO_TEST_CASE( testEmptyHumanReadableOption )
{
    const std::vector< std::string > options = {};
    const std::string optionString = echidna::caller::params::displayOptions( options );
    BOOST_CHECK_EQUAL( "", optionString );
}

BOOST_AUTO_TEST_CASE( testSingleHumanReadableStringOption )
{
    const std::vector< std::string > options = {"one"};
    const std::string optionString = echidna::caller::params::displayOptions( options );
    BOOST_CHECK_EQUAL( "one", optionString );
}

BOOST_AUTO_TEST_CASE( testTwoHumanReadableStringOptions )
{
    const std::vector< std::string > options = {"one", "two"};
    const std::string optionString = echidna::caller::params::displayOptions( options );
    BOOST_CHECK_EQUAL( "one or two", optionString );
}

BOOST_AUTO_TEST_CASE( testManyHumanReadableStringOptions )
{
    const std::vector< std::string > options = {"one", "two", "three"};
    const std::string optionString = echidna::caller::params::displayOptions( options );
    BOOST_CHECK_EQUAL( "one, two or three", optionString );
}

BOOST_AUTO_TEST_CASE( testSingleHumanReadableIntOption )
{
    const std::vector< int > options = {1};
    const std::string optionString = echidna::caller::params::displayOptions( options );
    BOOST_CHECK_EQUAL( "1", optionString );
}

BOOST_AUTO_TEST_CASE( testTwoHumanReadableIntOptions )
{
    const std::vector< int > options = {1, 2};
    const std::string optionString = echidna::caller::params::displayOptions( options );
    BOOST_CHECK_EQUAL( "1 or 2", optionString );
}

BOOST_AUTO_TEST_CASE( testManyHumanReadableIntOptions )
{
    const std::vector< int > options = {1, 2, 3};
    const std::string optionString = echidna::caller::params::displayOptions( options );
    BOOST_CHECK_EQUAL( "1, 2 or 3", optionString );
}
