#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/algorithm/string.hpp>

#include "vcf/filterDescription.hpp"
#include "utils/partition.hpp"

using output_t = std::vector< std::vector< int > >;

BOOST_AUTO_TEST_CASE( partitionShouldPartitionEmptyList )
{
    std::vector< int > input = {};
    output_t expectedOutput = {};
    output_t output = echidna::utils::functional::partition( input, []( int a, int b )
                                                             {
                                                                 return a == b;
                                                             } );

    BOOST_REQUIRE_EQUAL( output.size(), expectedOutput.size() );
    for ( std::size_t i = 0; i < output.size(); ++i )
    {
        BOOST_CHECK_EQUAL_COLLECTIONS( output[i].begin(), output[i].end(), expectedOutput[i].begin(),
                                       expectedOutput[i].end() );
    }
}

BOOST_AUTO_TEST_CASE( partitionShouldPartitionListWithOneElement )
{
    std::vector< int > input = {1};
    output_t expectedOutput = {{1}};
    output_t output = echidna::utils::functional::partition( input, []( int a, int b )
                                                             {
                                                                 return a == b;
                                                             } );

    BOOST_REQUIRE_EQUAL( output.size(), expectedOutput.size() );
    for ( std::size_t i = 0; i < output.size(); ++i )
    {
        BOOST_CHECK_EQUAL_COLLECTIONS( output[i].begin(), output[i].end(), expectedOutput[i].begin(),
                                       expectedOutput[i].end() );
    }
}

BOOST_AUTO_TEST_CASE( partitionShouldPreserveIntraListOrder )
{
    std::vector< int > input = {1, 2};
    output_t expectedOutput = {{1}, {2}};
    output_t output = echidna::utils::functional::partition( input, []( int a, int b )
                                                             {
                                                                 return a == b;
                                                             } );

    BOOST_REQUIRE_EQUAL( output.size(), expectedOutput.size() );
    for ( std::size_t i = 0; i < output.size(); ++i )
    {
        BOOST_CHECK_EQUAL_COLLECTIONS( output[i].begin(), output[i].end(), expectedOutput[i].begin(),
                                       expectedOutput[i].end() );
    }
}

BOOST_AUTO_TEST_CASE( partitionShouldPreserveOuterListOrder )
{
    std::vector< int > input = {1, 1};
    output_t expectedOutput = {{1, 1}};
    output_t output = echidna::utils::functional::partition( input, []( int a, int b )
                                                             {
                                                                 return a == b;
                                                             } );

    BOOST_REQUIRE_EQUAL( output.size(), expectedOutput.size() );
    for ( std::size_t i = 0; i < output.size(); ++i )
    {
        BOOST_CHECK_EQUAL_COLLECTIONS( output[i].begin(), output[i].end(), expectedOutput[i].begin(),
                                       expectedOutput[i].end() );
    }
}
