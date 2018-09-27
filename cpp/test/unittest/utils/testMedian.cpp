// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/median.hpp"

BOOST_AUTO_TEST_CASE( testMedianThrowsWithEmptyList )
{
    std::vector< double > input = {};
    BOOST_CHECK_THROW( echidna::utils::functional::median( input ), echidna::utils::echidna_exception );
}

BOOST_AUTO_TEST_CASE( testMedianWithOneElement )
{
    std::vector< double > input = {1.0};
    const auto expectedResult = 1.0;
    BOOST_CHECK_EQUAL( expectedResult, echidna::utils::functional::median( input ) );
}

BOOST_AUTO_TEST_CASE( testMedianWithTwoElements )
{
    std::vector< double > input = {3.0, 1.0};
    const auto expectedResult = 2.0;
    BOOST_CHECK_EQUAL( expectedResult, echidna::utils::functional::median( input ) );
}

BOOST_AUTO_TEST_CASE( testMedianWithTwoIntegerElements )
{
    std::vector< int > input = {2, 1};
    const auto expectedResult = 1;  // Mean of 2 & 1 is 1 for integers.
    BOOST_CHECK_EQUAL( expectedResult, echidna::utils::functional::median( input ) );
}

BOOST_AUTO_TEST_CASE( testMedianWithManyElements )
{
    std::vector< double > input = {8, 5, 9, 10, 1, 3, 4, 6, 7, 2};
    const auto expectedResult = 5.5;
    BOOST_CHECK_EQUAL( expectedResult, echidna::utils::functional::median( input ) );
}
