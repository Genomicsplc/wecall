#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/interval.hpp"

using echidna::utils::Interval;

BOOST_AUTO_TEST_CASE( emptyIntervalsShouldOverlapEachOther )
{
    BOOST_CHECK( Interval( 1, 1 ).overlaps( Interval( 1, 1 ) ) );
}

BOOST_AUTO_TEST_CASE( emptyIntervalDoesNotOverlapTouchingIntervalOnLeft )
{
    BOOST_CHECK( not Interval( 1, 2 ).overlaps( Interval( 2, 2 ) ) );
}

BOOST_AUTO_TEST_CASE( emptyIntervalDoesNotOverlapTouchingIntervalOnRight )
{
    BOOST_CHECK( not Interval( 1, 1 ).overlaps( Interval( 1, 2 ) ) );
}

BOOST_AUTO_TEST_CASE( shouldCombineNonoverlappingSortedIntervals )
{
    Interval intvl1( 1, 3 );
    Interval intvl2( 5, 10 );

    Interval expectedIntvl( 1, 10 );
    intvl1.combine( intvl2 );

    BOOST_CHECK_EQUAL( intvl1, expectedIntvl );
}

BOOST_AUTO_TEST_CASE( shouldCombineNonoverlappingNonsortedIntervals )
{
    Interval intvl1( 5, 10 );
    Interval intvl2( 1, 3 );

    Interval expectedIntvl( 1, 10 );
    intvl1.combine( intvl2 );

    BOOST_CHECK_EQUAL( intvl1, expectedIntvl );
}

BOOST_AUTO_TEST_CASE( shouldCombineOverlappingSortedIntervals )
{
    Interval intvl1( 1, 8 );
    Interval intvl2( 5, 10 );

    Interval expectedIntvl( 1, 10 );
    intvl1.combine( intvl2 );

    BOOST_CHECK_EQUAL( intvl1, expectedIntvl );
}
