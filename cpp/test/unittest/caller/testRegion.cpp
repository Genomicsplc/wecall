// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "caller/region.hpp"
#include "caller/regionUtils.hpp"
#include "utils/exceptions.hpp"
#include "io/fastaFile.hpp"

using wecall::utils::Interval;

using wecall::caller::Region;
using wecall::caller::SetRegions;

BOOST_AUTO_TEST_CASE( testGetPaddedByAmountProvided )
{
    BOOST_CHECK_EQUAL( Region( "1", 0, 10 ).getPadded( 5 ), Region( "1", -5, 15 ) );
    BOOST_CHECK_EQUAL( Region( "1", 0, 10 ).getPadded( -5 ), Region( "1", 5, 5 ) );
    BOOST_CHECK_THROW( Region( "1", 0, 11 ).getPadded( -6 ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( testEquality )
{
    Region region( "1", 0, 20 );
    BOOST_CHECK( region == region );
    BOOST_CHECK( Region( "1", 0, 10 ) == Region( "1", 0, 10 ) );

    BOOST_CHECK( not( region != region ) );
    BOOST_CHECK( not( Region( "1", 0, 10 ) != Region( "1", 0, 10 ) ) );

    BOOST_CHECK( not( Region( "1", 0, 10 ) == Region( "1", 5, 10 ) ) );
    BOOST_CHECK( not( Region( "1", 0, 10 ) == Region( "1", 0, 11 ) ) );
    BOOST_CHECK( not( Region( "1", 0, 10 ) == Region( "3", 0, 10 ) ) );

    BOOST_CHECK( Region( "1", 0, 10 ) != Region( "1", 5, 10 ) );
    BOOST_CHECK( Region( "1", 0, 10 ) != Region( "1", 0, 11 ) );
    BOOST_CHECK( Region( "1", 0, 10 ) != Region( "3", 0, 10 ) );
}

BOOST_AUTO_TEST_CASE( testContructorEquality )
{
    BOOST_CHECK( Region( "1", 0, 10 ) == Region( "1", Interval( 0, 10 ) ) );
}

BOOST_AUTO_TEST_CASE( testContains )
{
    BOOST_CHECK( Region( "1", 5, 10 ).contains( Region( "1", 7, 9 ) ) );
    BOOST_CHECK( Region( "1", 5, 10 ).contains( Region( "1", 5, 10 ) ) );
    BOOST_CHECK( not Region( "1", 5, 10 ).contains( Region( "2", 5, 10 ) ) );
    BOOST_CHECK( not Region( "1", 5, 10 ).contains( Region( "1", 4, 10 ) ) );
    BOOST_CHECK( not Region( "1", 5, 10 ).contains( Region( "1", 5, 11 ) ) );
    BOOST_CHECK( Region( "1", 5, 10 ).contains( "1", 5 ) );
    BOOST_CHECK( Region( "1", 5, 10 ).contains( "1", 9 ) );
    BOOST_CHECK( not Region( "1", 5, 10 ).contains( "1", 10 ) );
    BOOST_CHECK( not Region( "1", 5, 10 ).contains( "1", 4 ) );
    BOOST_CHECK( not Region( "1", 5, 10 ).contains( "2", 6 ) );
}

BOOST_AUTO_TEST_CASE( shouldOverlapEmptyRegion )
{
    BOOST_CHECK( Region( "", 1, 1 ).overlaps( Region( "", 1, 1 ) ) );
    BOOST_CHECK( Region( "", 1, 1 ).overlaps( Region( "", 1, 1 ) ) );

    BOOST_CHECK( Region( "", 0, 2 ).overlaps( Region( "", 1, 1 ) ) );
    BOOST_CHECK( Region( "", 1, 1 ).overlaps( Region( "", 0, 2 ) ) );

    BOOST_CHECK( not Region( "", 0, 1 ).overlaps( Region( "", 1, 1 ) ) );
    BOOST_CHECK( not Region( "", 1, 1 ).overlaps( Region( "", 0, 1 ) ) );

    BOOST_CHECK( not Region( "", 1, 2 ).overlaps( Region( "", 1, 1 ) ) );
    BOOST_CHECK( not Region( "", 1, 1 ).overlaps( Region( "", 1, 2 ) ) );
}

BOOST_AUTO_TEST_CASE( shouldThrowWhenCombiningRegionsOnDifferentContigs )
{
    Region regionChr1( "1", 5, 10 );
    Region regionChr2( "2", 8, 15 );
    BOOST_CHECK_THROW( regionChr1.combine( regionChr2 ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( shouldCombineNonOverlappingRegionsWithSameContig )
{
    Region regionChr1( "1", 5, 10 );
    Region regionChr2( "1", 15, 20 );
    regionChr1.combine( regionChr2 );
    BOOST_CHECK_EQUAL( regionChr1, Region( "1", 5, 20 ) );
}

BOOST_AUTO_TEST_CASE( shouldOutputStringInSamtoolsFormat )
{
    Region region( "2", 0, 10 );
    std::string expectedString = "2:0-10";
    BOOST_CHECK_EQUAL( expectedString, region.toString() );
}

BOOST_AUTO_TEST_CASE( testRegionSetFillEmptySet )
{
    SetRegions setRegions;

    BOOST_CHECK_EQUAL( setRegions.size(), 0 );

    setRegions.fill( 100 );
    BOOST_CHECK_EQUAL( setRegions.size(), 0 );
}

BOOST_AUTO_TEST_CASE( testRegionSetFillOneItem )
{
    SetRegions setRegions;
    setRegions.insert( Region( "1", 2, 3 ) );

    BOOST_CHECK_EQUAL( setRegions.size(), 1 );

    setRegions.fill( 100 );
    BOOST_CHECK_EQUAL( setRegions.size(), 1 );
    BOOST_CHECK_EQUAL( setRegions.getSpan(), Region( "1", 2, 3 ) );
}

BOOST_AUTO_TEST_CASE( testRegionSetFill )
{
    SetRegions setRegions;
    setRegions.insert( Region( "1", 0, 0 ) );
    setRegions.insert( Region( "1", 2, 3 ) );

    BOOST_CHECK_EQUAL( setRegions.size(), 2 );

    setRegions.fill( 1 );
    BOOST_CHECK_EQUAL( setRegions.size(), 2 );
    BOOST_CHECK_EQUAL( setRegions.getSpan(), Region( "1", 0, 3 ) );

    setRegions.fill( 2 );
    BOOST_CHECK_EQUAL( setRegions.size(), 1 );
    BOOST_CHECK_EQUAL( setRegions.getSpan(), Region( "1", 0, 3 ) );
}
