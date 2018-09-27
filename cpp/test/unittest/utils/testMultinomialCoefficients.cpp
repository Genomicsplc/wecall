// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "utils/multinomialCoefficients.hpp"

// n=1; r=1

BOOST_AUTO_TEST_CASE( multinomialCoefficient_1 )
{
    std::vector< unsigned int > input = {1};
    BOOST_CHECK_EQUAL( 1, multinomial_coefficient( input ) );
}

// n=2; r=2

BOOST_AUTO_TEST_CASE( multinomialCoefficient_0_2 )
{
    std::vector< unsigned int > input = {0, 2};
    BOOST_CHECK_EQUAL( 1, multinomial_coefficient( input ) );
}

BOOST_AUTO_TEST_CASE( multinomialCoefficient_1_1 )
{
    std::vector< unsigned int > input = {1, 1};
    BOOST_CHECK_EQUAL( 2, multinomial_coefficient( input ) );
}

BOOST_AUTO_TEST_CASE( multinomialCoefficient_2_0 )
{
    std::vector< unsigned int > input = {2, 0};
    BOOST_CHECK_EQUAL( 1, multinomial_coefficient( input ) );
}

// n=3; r=3 (partial)

BOOST_AUTO_TEST_CASE( multinomialCoefficient_3_0_0 )
{
    std::vector< unsigned int > input = {3, 0, 0};
    BOOST_CHECK_EQUAL( 1, multinomial_coefficient( input ) );
}

BOOST_AUTO_TEST_CASE( multinomialCoefficient_2_1_0 )
{
    std::vector< unsigned int > input = {2, 1, 0};
    BOOST_CHECK_EQUAL( 3, multinomial_coefficient( input ) );
}

BOOST_AUTO_TEST_CASE( multinomialCoefficient_2_0_1 )
{
    std::vector< unsigned int > input = {2, 0, 1};
    BOOST_CHECK_EQUAL( 3, multinomial_coefficient( input ) );
}

BOOST_AUTO_TEST_CASE( multinomialCoefficient_1_1_1 )
{
    std::vector< unsigned int > input = {1, 1, 1};
    BOOST_CHECK_EQUAL( 6, multinomial_coefficient( input ) );
}

BOOST_AUTO_TEST_CASE( multinomialCoefficient_0_2_1 )
{
    std::vector< unsigned int > input = {0, 2, 1};
    BOOST_CHECK_EQUAL( 3, multinomial_coefficient( input ) );
}

BOOST_AUTO_TEST_CASE( multinomialCoefficient_0_1_2 )
{
    std::vector< unsigned int > input = {0, 1, 2};
    BOOST_CHECK_EQUAL( 3, multinomial_coefficient( input ) );
}

BOOST_AUTO_TEST_CASE( multinomialCoefficient_0_0_3 )
{
    std::vector< unsigned int > input = {0, 0, 3};
    BOOST_CHECK_EQUAL( 1, multinomial_coefficient( input ) );
}
