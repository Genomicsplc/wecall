// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "caller/region.hpp"

#include "utils/referenceSequence.hpp"

using namespace wecall;
using wecall::utils::ReferenceSequence;
using wecall::utils::BasePairSequence;
using wecall::utils::Interval;
using wecall::caller::Region;

BOOST_AUTO_TEST_CASE( testConstructorSanity )
{
    caller::Region r( "1", utils::Interval( 0, 1 ) );

    BOOST_CHECK_THROW( utils::ReferenceSequence( r, "hoho" ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( testConstructorInsanity )
{
    caller::Region r( "1", utils::Interval( 0, 1 ) );

    BOOST_CHECK_EQUAL( utils::ReferenceSequence( r, "a" ).region(), r );
}

BOOST_AUTO_TEST_CASE( testSequenceIsGettableOutable )
{
    caller::Region r( "1", utils::Interval( 0, 1 ) );
    const wecall::utils::BasePairSequence seq( "a" );
    BOOST_CHECK_EQUAL( utils::ReferenceSequence( r, seq ).sequence(), seq );
}

BOOST_AUTO_TEST_CASE( testSubseq )
{
    caller::Region r( "1", utils::Interval( 101, 103 ) );
    caller::Region subRegion( "1", utils::Interval( 101, 102 ) );
    const wecall::utils::BasePairSequence seq( "ab" );
    const wecall::utils::BasePairSequence expectedSeq( "a" );
    BOOST_CHECK_EQUAL( utils::ReferenceSequence( r, seq ).subseq( subRegion ).region(), subRegion );
    BOOST_CHECK_EQUAL( utils::ReferenceSequence( r, seq ).subseq( subRegion ).sequence(), expectedSeq );
}

BOOST_AUTO_TEST_CASE( testSubseqRange )
{
    caller::Region r( "1", utils::Interval( 101, 103 ) );
    caller::Region subRegion( "1", utils::Interval( 101, 102 ) );
    const wecall::utils::BasePairSequence seq( "ab" );
    const wecall::utils::BasePairSequence expectedSeq( "a" );
    auto range = utils::ReferenceSequence( r, seq ).subseqRange( subRegion );
    BOOST_CHECK_EQUAL( std::string( range.first, range.second ), expectedSeq.str() );
}

BOOST_AUTO_TEST_CASE( testAccessors )
{
    std::string const contig = "1";
    caller::Region r( contig, utils::Interval( 100, 104 ) );
    wecall::utils::ReferenceSequence s( r, "ACTG" );

    BOOST_CHECK_EQUAL( s.size(), 4 );
    BOOST_CHECK_EQUAL( s.start(), 100 );
    BOOST_CHECK_EQUAL( s.end(), 104 );
    BOOST_CHECK_EQUAL( s.contig(), contig );
}

BOOST_AUTO_TEST_CASE( testGetPaddedReference )
{
    caller::Region r( "1", utils::Interval( 10, 11 ) );
    const wecall::utils::BasePairSequence seq( "a" );

    BOOST_CHECK_EQUAL(
        utils::ReferenceSequence( r, "a" ).getPadded( caller::Region( "1", utils::Interval( 8, 14 ) ) ).region(),
        caller::Region( "1", utils::Interval( 8, 14 ) ) );
    BOOST_CHECK_EQUAL(
        utils::ReferenceSequence( r, "a" ).getPadded( caller::Region( "1", utils::Interval( 8, 14 ) ) ).sequence(),
        "NNaNNN" );
}

BOOST_AUTO_TEST_CASE( testBasicForwardIterators )
{
    const ReferenceSequence referenceSequence( Region( "1", 1, 11 ), "0123456789" );

    auto it = referenceSequence.cbegin();
    auto itEnd = referenceSequence.cend();
    BOOST_REQUIRE_EQUAL( std::distance( it, itEnd ), 10 );

    {
        std::string obtainedSequenceStr( 10, ' ' );
        std::size_t i = 0;
        for ( ; it != itEnd; ++it )
        {
            obtainedSequenceStr[i] = *it;
            ++i;
        }
        wecall::utils::BasePairSequence obtainedSequence = obtainedSequenceStr;

        BOOST_CHECK_EQUAL( obtainedSequence, referenceSequence.sequence() );
    }

    BOOST_CHECK_THROW( referenceSequence.getRangeForwardIterators( Interval( 0, 5 ) ),
                       wecall::utils::wecall_exception );
    BOOST_CHECK_THROW( referenceSequence.getRangeForwardIterators( Interval( 5, 12 ) ),
                       wecall::utils::wecall_exception );

    BOOST_CHECK_EQUAL( *referenceSequence.getRangeForwardIterators( Interval( 1, 5 ) ).first, '0' );
    BOOST_CHECK_EQUAL( *referenceSequence.getRangeForwardIterators( Interval( 1, 5 ) ).second, '4' );
    BOOST_CHECK_EQUAL( *referenceSequence.getRangeForwardIterators( Interval( 2, 6 ) ).first, '1' );
    BOOST_CHECK_EQUAL( *referenceSequence.getRangeForwardIterators( Interval( 2, 6 ) ).second, '5' );
    BOOST_CHECK_EQUAL( *referenceSequence.getRangeForwardIterators( Interval( 10, 11 ) ).first, '9' );
    BOOST_CHECK( referenceSequence.getRangeForwardIterators( Interval( 10, 11 ) ).second == referenceSequence.cend() );
    BOOST_CHECK( referenceSequence.getRangeForwardIterators( Interval( 11, 11 ) ).first == referenceSequence.cend() );
    BOOST_CHECK( referenceSequence.getRangeForwardIterators( Interval( 11, 11 ) ).second == referenceSequence.cend() );

    const auto range = referenceSequence.getRangeForwardIterators( Interval( 3, 7 ) );
    BOOST_CHECK_EQUAL( std::distance( range.first, range.second ), 4 );
    {
        std::string obtainedSequenceStr( 4, ' ' );
        std::size_t i = 0;
        it = range.first;
        itEnd = range.second;
        for ( ; it != itEnd; ++it )
        {
            obtainedSequenceStr[i] = *it;
            ++i;
        }
        wecall::utils::BasePairSequence obtainedSequence = obtainedSequenceStr;
        BOOST_CHECK_EQUAL( obtainedSequenceStr, "2345" );
    }
}

BOOST_AUTO_TEST_CASE( testBasicReverseIterators )
{
    const ReferenceSequence referenceSequence( Region( "1", 1, 11 ), "0123456789" );

    auto it = referenceSequence.crbegin();
    auto itEnd = referenceSequence.crend();
    BOOST_REQUIRE_EQUAL( std::distance( it, itEnd ), 10 );

    {
        std::string obtainedSequenceStr( 10, ' ' );
        std::size_t i = 0;
        for ( ; it != itEnd; ++it )
        {
            obtainedSequenceStr[i] = *it;
            ++i;
        }
        wecall::utils::BasePairSequence obtainedSequence = obtainedSequenceStr;
        BOOST_CHECK_EQUAL( obtainedSequence, "9876543210" );
    }

    BOOST_CHECK_THROW( referenceSequence.getRangeReverseIterators( Interval( 0, 5 ) ),
                       wecall::utils::wecall_exception );
    BOOST_CHECK_THROW( referenceSequence.getRangeReverseIterators( Interval( 5, 12 ) ),
                       wecall::utils::wecall_exception );

    BOOST_CHECK( referenceSequence.getRangeReverseIterators( Interval( 1, 5 ) ).second == referenceSequence.crend() );
    BOOST_CHECK_EQUAL( *referenceSequence.getRangeReverseIterators( Interval( 1, 5 ) ).first, '3' );

    BOOST_CHECK_EQUAL( *referenceSequence.getRangeReverseIterators( Interval( 2, 6 ) ).first, '4' );
    BOOST_CHECK_EQUAL( *referenceSequence.getRangeReverseIterators( Interval( 2, 6 ) ).second, '0' );

    BOOST_CHECK_EQUAL( *referenceSequence.getRangeReverseIterators( Interval( 10, 11 ) ).first, '9' );
    BOOST_CHECK_EQUAL( *referenceSequence.getRangeReverseIterators( Interval( 10, 11 ) ).second, '8' );

    BOOST_CHECK_EQUAL( *referenceSequence.getRangeReverseIterators( Interval( 11, 11 ) ).first, '9' );
    BOOST_CHECK_EQUAL( *referenceSequence.getRangeReverseIterators( Interval( 11, 11 ) ).second, '9' );

    const auto range = referenceSequence.getRangeReverseIterators( Interval( 3, 7 ) );
    BOOST_CHECK_EQUAL( std::distance( range.first, range.second ), 4 );
    {
        std::string obtainedSequenceStr( 4, ' ' );
        std::size_t i = 0;
        it = range.first;
        itEnd = range.second;
        for ( ; it != itEnd; ++it )
        {
            obtainedSequenceStr[i] = *it;
            ++i;
        }
        wecall::utils::BasePairSequence obtainedSequence = obtainedSequenceStr;
        BOOST_CHECK_EQUAL( obtainedSequenceStr, "5432" );
    }
}
