// All content Copyright (C) 2018 Genomics plc
#include "alignment/mmHelpers.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <iostream>

BOOST_AUTO_TEST_CASE( testMmRepeatContructor )
{
    auto s = short( 0x12345678 );
    const short_array8 a( s );
    for ( auto i = 0; i < 8; i++ )
        BOOST_CHECK_EQUAL( a[i], s );
}

#ifndef WECALL_MM_HELPERS2_H
BOOST_AUTO_TEST_CASE( testMmContructorFromMmSet )
{
    __m128i s = _mm_set_epi16( 0, 1, 2, 3, 4, 5, 6, 7 );
    const short_array8 a( s );
    for ( auto i = 0; i < 8; i++ )
        BOOST_CHECK_EQUAL( a[i], 7 - i );
}
#endif

BOOST_AUTO_TEST_CASE( testMmContructorFromShorts )
{
    const short_array8 a( 0, 1, 2, 3, 4, 5, 6, 7 );
    for ( auto i = 0; i < 8; i++ )
        BOOST_CHECK_EQUAL( a[i], i );
}

BOOST_AUTO_TEST_CASE( testMmDefaultContructor )
{
    const short_array8 a;
    for ( auto i = 0; i < 8; i++ )
        BOOST_CHECK_EQUAL( a[i], 0 );
}

BOOST_AUTO_TEST_CASE( testMmCopyContructor )
{
    const short_array8 a1( 0, 1, 2, 3, 4, 5, 6, 7 );
    const short_array8 a2( a1 );
    for ( auto i = 0; i < 8; i++ )
        BOOST_CHECK_EQUAL( a1[i], a2[i] );
}

BOOST_AUTO_TEST_CASE( testMmContructWithLambda )
{
    short mult = 2;
    const short_array8 a1( [mult]( short i )
                           {
                               return i * mult;
                           } );
    for ( auto i = 0; i < 8; i++ )
        BOOST_CHECK_EQUAL( a1[i], i * mult );
}

BOOST_AUTO_TEST_CASE( testMmNonConstReferenceRead )
{
    short_array8 a( 0, 1, 2, 3, 4, 5, 6, 7 );
    for ( auto i = 0; i < 8; i++ )
        BOOST_CHECK_EQUAL( a[i], i );
}

BOOST_AUTO_TEST_CASE( testMmNonConstReferenceSet )
{
    short_array8 a( short( 0xffff ) );
    auto & ca = a;

    a[0] = 0;
    a[1] = 1;
    a[2] = 2;
    a[3] = 3;
    a[4] = 4;
    a[5] = 5;
    a[6] = 6;
    a[7] = 7;

    for ( short i = 0; i < 8; i++ )
        BOOST_CHECK_EQUAL( ca[i], i );
    for ( auto i = 0; i < 8; i++ )
        BOOST_CHECK_EQUAL( a[i], i );
}

BOOST_AUTO_TEST_CASE( testMmAssign )
{
    const short_array8 a1( 0, 1, 2, 3, 4, 5, 6, 7 );
    short_array8 a2( short( 0xffff ) );
    a2 = a1;
    for ( auto i = 0; i < 8; i++ )
        BOOST_CHECK_EQUAL( a1[i], a2[i] );
}

BOOST_AUTO_TEST_CASE( testMmMin )
{
    const short_array8 a1( 1, 2, 99, -4, 0, 6, 7, 8 );
    const short_array8 a2( 1, 3, 2, -4, 5, 0, 7, 8 );
    const short_array8 ex( 1, 2, 2, -4, 0, 0, 7, 8 );
    const short_array8 m = min( a1, a2 );
    for ( auto i = 0; i < 8; i++ )
        BOOST_CHECK_EQUAL( m[i], ex[i] );
}

BOOST_AUTO_TEST_CASE( testMmAndNot )
{
    const short_array8 a1( short( 0xffff ), 0x0000, short( 0xff00 ), 1, 1, 0, 0, 0 );
    const short_array8 a2( short( 0xffff ), short( 0xffff ), 0x00ff, 1, 0, 1, 0, 0 );
    const short_array8 ex( 0x0000, short( 0xffff ), 0x00ff, 0, 0, 1, 0, 0 );
    const short_array8 an = andnot( a1, a2 );
    for ( auto i = 0; i < 8; i++ )
        BOOST_CHECK_EQUAL( an[i], ex[i] );
}

BOOST_AUTO_TEST_CASE( testMmCmpEq )
{
    const short f = short( 0xffff );
    const short_array8 a1( 01, 00, 01, f, 0x0000, 0, 0, 0 );
    const short_array8 a2( 01, 01, 00, f, 0x0000, 1, 0, 0 );
    const short_array8 ex( -1, 00, 00, f, f, 0, f, f );
    const short_array8 ac = cmpeq( a1, a2 );
    for ( auto i = 0; i < 8; i++ )
    {
        BOOST_CHECK_EQUAL( ac[i], ex[i] );
    }
}

BOOST_AUTO_TEST_CASE( testMmPlusOperator )
{
    const short_array8 a1( 1, short( 0xffff ), short( 0xff00 ), 0x0001, 0x0000, 0, -5, 3 );
    const short_array8 a2( 2, short( 0xffff ), 0x00ff, 0x0001, 0x0000, 1, 5, -99 );
    const short_array8 ex( 3, short( 0xffff + 0xffff ), short( 0xffff ), 0x2, 0x0, 0x0001, 0, -96 );
    const short_array8 ac = a1 + a2;
    for ( auto i = 0; i < 8; i++ )
        BOOST_CHECK_EQUAL( ac[i], ex[i] );
}

BOOST_AUTO_TEST_CASE( testMmBitwiseAnd )
{
    const short_array8 a1( short( 0xffff ), 0x0000, short( 0xff00 ), 0x0001, 0x0000, 0, 0, 3 );
    const short_array8 a2( short( 0xffff ), short( 0xffff ), 0x00ff, 0x0001, 0x0000, 1, 0, 6 );
    const short_array8 ex( short( 0xffff ), 0x0000, 0x0000, 0x0001, 0x0000, 0, 0, 2 );
    const short_array8 ac = a1 & a2;
    for ( auto i = 0; i < 8; i++ )
        BOOST_CHECK_EQUAL( ac[i], ex[i] );
}

BOOST_AUTO_TEST_CASE( testMmBitwiseOr )
{
    const short_array8 a1( short( 0xffff ), 0x0000, short( 0xff00 ), 0x0001, 0x0000, 0, 0, 3 );
    const short_array8 a2( short( 0xffff ), short( 0xffff ), 0x00ff, 0x0001, 0x0000, 1, 0, 6 );
    const short_array8 ex( short( 0xffff ), short( 0xffff ), short( 0xffff ), 0x0001, 0x0000, 1, 0, 7 );
    const short_array8 ac = a1 | a2;
    for ( auto i = 0; i < 8; i++ )
        BOOST_CHECK_EQUAL( ac[i], ex[i] );
}

BOOST_AUTO_TEST_CASE( testMmBitwiseShiftLeft )
{
    const short_array8 a1( short( 0xffff ), 0x0000, short( 0xff00 ), 0x0001, 0x0000, 0, 0, 3 );
    const short_array8 ex( short( 0xfffc ), 0x0000, short( 0xfc00 ), 0x0004, 0x0000, 0, 0, 0xc );

    const short_array8 ac = a1 << 2;
    for ( auto i = 0; i < 8; i++ )
        BOOST_CHECK_EQUAL( ac[i], ex[i] );
}

BOOST_AUTO_TEST_CASE( testMmBitwiseShiftRight )
{
    const short_array8 a1( short( 0xffff ), 0x0000, short( 0xff00 ), 0x0001, 0x0000, 0, 0, 0x8 );
    const short_array8 ex( 0x3fff, 0x0000, 0x3fc0, 0x0000, 0x0000, 0, 0, 0x2 );

    const short_array8 ac = a1 >> 2;
    for ( auto i = 0; i < 8; i++ )
    {
        BOOST_CHECK_EQUAL( ac[i], ex[i] );
    }
}
