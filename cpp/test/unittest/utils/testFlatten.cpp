// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/algorithm/string.hpp>

#include "vcf/filterDescription.hpp"
#include "utils/flatten.hpp"

BOOST_AUTO_TEST_CASE( shouldFlattenEmptyList )
{
    std::vector< std::vector< int > > input = {{}};
    std::vector< int > expectedOutput = {};
    std::vector< int > output = echidna::utils::functional::flatten( input );
    BOOST_CHECK_EQUAL_COLLECTIONS( output.begin(), output.end(), expectedOutput.begin(), expectedOutput.end() );
}

BOOST_AUTO_TEST_CASE( shouldFlattenEmptyInnerList )
{
    std::vector< std::vector< int > > input = {};
    std::vector< int > expectedOutput = {};
    std::vector< int > output = echidna::utils::functional::flatten( input );
    BOOST_CHECK_EQUAL_COLLECTIONS( output.begin(), output.end(), expectedOutput.begin(), expectedOutput.end() );
}

BOOST_AUTO_TEST_CASE( shouldFlattenListWithOneElement )
{
    std::vector< std::vector< int > > input = {{1}};
    std::vector< int > expectedOutput = {1};
    std::vector< int > output = echidna::utils::functional::flatten( input );
    BOOST_CHECK_EQUAL_COLLECTIONS( output.begin(), output.end(), expectedOutput.begin(), expectedOutput.end() );
}

BOOST_AUTO_TEST_CASE( shouldPreserveIntraListOrder )
{
    std::vector< std::vector< int > > input = {{1, 2}};
    std::vector< int > expectedOutput = {1, 2};
    std::vector< int > output = echidna::utils::functional::flatten( input );
    BOOST_CHECK_EQUAL_COLLECTIONS( output.begin(), output.end(), expectedOutput.begin(), expectedOutput.end() );
}

BOOST_AUTO_TEST_CASE( shouldPreserveOuterListOrder )
{
    std::vector< std::vector< int > > input = {{1}, {2}};
    std::vector< int > expectedOutput = {1, 2};
    std::vector< int > output = echidna::utils::functional::flatten( input );
    BOOST_CHECK_EQUAL_COLLECTIONS( output.begin(), output.end(), expectedOutput.begin(), expectedOutput.end() );
}
