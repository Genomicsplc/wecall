#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/sequence.hpp"

using namespace echidna;

BOOST_AUTO_TEST_CASE( testBasePairSequenceConstructors )
{
    utils::BasePairSequence a( "A" );
    utils::BasePairSequence b = a;
    utils::BasePairSequence c( std::string( "A" ) );
    utils::BasePairSequence d( 1, 'A' );

    BOOST_CHECK_EQUAL( a, a );
    BOOST_CHECK_EQUAL( a, b );
    BOOST_CHECK_EQUAL( a, c );
    BOOST_CHECK_EQUAL( a, d );
}

BOOST_AUTO_TEST_CASE( testBasePairSequenceEquality )
{
    BOOST_CHECK( utils::BasePairSequence( "HOHO" ) == utils::BasePairSequence( "HOHO" ) );
    BOOST_CHECK( utils::BasePairSequence( "HEHE" ) != utils::BasePairSequence( "HOHO" ) );
}

BOOST_AUTO_TEST_CASE( testBasePairSequenceComparison )
{
    BOOST_CHECK( utils::BasePairSequence( "HEHE" ) < utils::BasePairSequence( "HOHO" ) );
    BOOST_CHECK( utils::BasePairSequence( "HEHE" ) <= utils::BasePairSequence( "HOHO" ) );

    BOOST_CHECK( utils::BasePairSequence( "HOHO" ) > utils::BasePairSequence( "HEHE" ) );
    BOOST_CHECK( utils::BasePairSequence( "HOHO" ) >= utils::BasePairSequence( "HEHE" ) );
}

BOOST_AUTO_TEST_CASE( testBasePairAccessors )
{
    utils::BasePairSequence seq( "ATCG" );
    BOOST_CHECK_EQUAL( seq[0], 'A' );
    BOOST_CHECK_EQUAL( seq.at( 0 ), 'A' );

    BOOST_CHECK_EQUAL( seq[1], 'T' );

    BOOST_CHECK_THROW( seq.at( 4 ), std::exception );
}

BOOST_AUTO_TEST_CASE( testBasePairConversionToString )
{
    utils::BasePairSequence seq( "ATCG" );
    BOOST_CHECK_EQUAL( seq.str(), std::string( "ATCG" ) );

    BOOST_CHECK_EQUAL( seq.size(), 4 );
}

BOOST_AUTO_TEST_CASE( testBasePairSequenceSubstr )
{
    utils::BasePairSequence seq( "ATCG" );
    utils::BasePairSequence a( "A" ), t( "T" ), c( "C" ), g( "G" );
    BOOST_CHECK_EQUAL( seq.substr( 0, 1 ), a );
    BOOST_CHECK_EQUAL( seq.substr( 1, 1 ), t );
    BOOST_CHECK_EQUAL( seq.substr( 2, 1 ), c );
    BOOST_CHECK_EQUAL( seq.substr( 3, 1 ), g );

    BOOST_CHECK_EQUAL( seq.substr( 4, 1 ), utils::BasePairSequence() );

    BOOST_CHECK_THROW( seq.substr( -1, 1 ), echidna::utils::echidna_exception );
}

BOOST_AUTO_TEST_CASE( testBasePairSequenceContains )
{
    utils::BasePairSequence seq( "ATCG" );
    BOOST_CHECK( not seq.contains( 'N' ) );
    BOOST_CHECK( seq.contains( 'A' ) );
    BOOST_CHECK( seq.contains( 'C' ) );
    BOOST_CHECK( seq.contains( 'T' ) );
    BOOST_CHECK( seq.contains( 'G' ) );
}

BOOST_AUTO_TEST_CASE( testBasePairSequenceForwardIterator )
{
    utils::BasePairSequence a( "AT" );

    auto it = a.cbegin();
    BOOST_CHECK_EQUAL( *it, 'A' );
    ++it;
    BOOST_CHECK_EQUAL( *it, 'T' );
    ++it;
    BOOST_CHECK( it == a.cend() );
}

BOOST_AUTO_TEST_CASE( testBasePairSequenceReverseIterator )
{
    utils::BasePairSequence a( "AT" );

    auto it = a.crbegin();
    BOOST_CHECK_EQUAL( *it, 'T' );
    ++it;
    BOOST_CHECK_EQUAL( *it, 'A' );
    ++it;
    BOOST_CHECK( it == a.crend() );
}

BOOST_AUTO_TEST_CASE( testBasePairSequenceReverseIteratorOnSubStr )
{
    utils::BasePairSequence a( "1234" );
    auto b = a.substr( 1, 2 );

    auto it = b.crbegin();
    BOOST_CHECK_EQUAL( *it, '3' );
    ++it;
    BOOST_CHECK_EQUAL( *it, '2' );
    ++it;
    BOOST_CHECK( it == b.crend() );
}

BOOST_AUTO_TEST_CASE( testStreamOperator )
{
    std::string const seqStr = "ATCG";
    utils::BasePairSequence seq( seqStr );

    std::stringstream out;

    out << seq;

    BOOST_CHECK_EQUAL( out.str(), seqStr );
}

BOOST_AUTO_TEST_CASE( testStringAddition )
{
    utils::BasePairSequence A( "G" );
    utils::BasePairSequence C( "T" );

    BOOST_CHECK_EQUAL( ( A + C ).str(), "GT" );
}

BOOST_AUTO_TEST_CASE( testLeftTrimmedSequenceNoGapChars )
{
    utils::BasePairSequence A( "GTC" );
    BOOST_CHECK_EQUAL( A.leftTrimmed(), A );
}

BOOST_AUTO_TEST_CASE( testRightTrimmedSequenceNoGapChars )
{
    utils::BasePairSequence A( "GTC" );
    BOOST_CHECK_EQUAL( A.rightTrimmed(), A );
}

BOOST_AUTO_TEST_CASE( testLeftTrimmedSequence )
{
    utils::BasePairSequence A( "NNNGTCNNN" );
    BOOST_CHECK_EQUAL( A.leftTrimmed().str(), "GTCNNN" );
}

BOOST_AUTO_TEST_CASE( testRightTrimmedSequence )
{
    utils::BasePairSequence A( "NNNGTCNNN" );
    BOOST_CHECK_EQUAL( A.rightTrimmed().str(), "NNNGTC" );
}

BOOST_AUTO_TEST_CASE( testDefaultContructor )
{
    utils::BasePairSequence A;
    utils::BasePairSequence B;
    BOOST_CHECK_EQUAL( A.str(), std::string() );
    BOOST_CHECK_EQUAL( A.size(), 0 );
    BOOST_CHECK( A.cbegin() == A.cend() );
    BOOST_CHECK( A.crbegin() == A.crend() );
    BOOST_CHECK_EQUAL( A, B );
}

BOOST_AUTO_TEST_CASE( testRepeatContructor )
{
    utils::BasePairSequence A( 10, 'A' );
    BOOST_CHECK_EQUAL( A.str(), std::string( 10, 'A' ) );
}