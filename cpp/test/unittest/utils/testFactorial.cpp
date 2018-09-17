#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "utils/multinomialCoefficients.hpp"

BOOST_AUTO_TEST_CASE( factorial_0 ) { BOOST_CHECK_EQUAL( 1, factorial( 0 ) ); }

BOOST_AUTO_TEST_CASE( factorial_1 ) { BOOST_CHECK_EQUAL( 1, factorial( 1 ) ); }

BOOST_AUTO_TEST_CASE( factorial_2 ) { BOOST_CHECK_EQUAL( 2, factorial( 2 ) ); }

BOOST_AUTO_TEST_CASE( factorial_3 ) { BOOST_CHECK_EQUAL( 6, factorial( 3 ) ); }
