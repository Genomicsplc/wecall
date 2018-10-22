// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "variant/type/breakpoint.hpp"
#include "utils/exceptions.hpp"

#include "caller/region.hpp"

using wecall::utils::BasePairSequence;
using wecall::variant::Breakpoint;
using wecall::variant::Variant;
using wecall::variant::BreakpointLocus;
using wecall::variant::breakpointPtrComp;
using wecall::variant::breakpointLocusPtrComp;
using wecall::utils::ReferenceSequence;
using wecall::caller::Region;

BOOST_AUTO_TEST_CASE( testConstructionOfBreakpoint )
{
    const Breakpoint startBreakpoint( "1", 20, true, "ATCG" );
    BOOST_CHECK_EQUAL( startBreakpoint.contig(), "1" );
    BOOST_CHECK_EQUAL( startBreakpoint.sequence(), "ATCG" );
    BOOST_CHECK_EQUAL( startBreakpoint.isStartBreakpoint(), true );

    const Breakpoint endBreakpoint( "1", 20, false, "ATCG" );
    BOOST_CHECK_EQUAL( endBreakpoint.contig(), "1" );
    BOOST_CHECK_EQUAL( endBreakpoint.sequence(), "ATCG" );
    BOOST_CHECK_EQUAL( endBreakpoint.isStartBreakpoint(), false );

    BOOST_CHECK( not( startBreakpoint == endBreakpoint ) );
}

BOOST_AUTO_TEST_CASE( testSortingOfBreakpoints )
{
    BOOST_CHECK( breakpointPtrComp()( std::make_shared< Breakpoint >( "a", 20, true, "A" ),
                                      std::make_shared< Breakpoint >( "b", 20, true, "A" ) ) );
    BOOST_CHECK( breakpointPtrComp()( std::make_shared< Breakpoint >( "1", 20, true, "A" ),
                                      std::make_shared< Breakpoint >( "1", 21, true, "A" ) ) );
    BOOST_CHECK( breakpointPtrComp()( std::make_shared< Breakpoint >( "1", 20, false, "A" ),
                                      std::make_shared< Breakpoint >( "1", 20, true, "A" ) ) );
    BOOST_CHECK( breakpointPtrComp()( std::make_shared< Breakpoint >( "1", 20, true, "B" ),
                                      std::make_shared< Breakpoint >( "1", 20, true, "AA" ) ) );
    BOOST_CHECK( breakpointPtrComp()( std::make_shared< Breakpoint >( "1", 20, true, "AA" ),
                                      std::make_shared< Breakpoint >( "1", 20, true, "AB" ) ) );
}

BOOST_AUTO_TEST_CASE( testSortingOfBreakpointLoci )
{
    // Sorts by contig
    BOOST_CHECK( breakpointLocusPtrComp()( std::make_shared< BreakpointLocus >( "a", 20, true ),
                                           std::make_shared< BreakpointLocus >( "b", 20, true ) ) );

    BOOST_CHECK( not breakpointLocusPtrComp()( std::make_shared< BreakpointLocus >( "b", 20, true ),
                                               std::make_shared< BreakpointLocus >( "a", 20, true ) ) );
    // Then by position
    BOOST_CHECK( breakpointLocusPtrComp()( std::make_shared< BreakpointLocus >( "1", 20, true ),
                                           std::make_shared< BreakpointLocus >( "1", 21, true ) ) );

    BOOST_CHECK( not breakpointLocusPtrComp()( std::make_shared< BreakpointLocus >( "1", 21, true ),
                                               std::make_shared< BreakpointLocus >( "1", 20, true ) ) );

    // Then putting 'end' loci before 'start' loci.
    BOOST_CHECK( breakpointLocusPtrComp()( std::make_shared< BreakpointLocus >( "1", 20, false ),
                                           std::make_shared< BreakpointLocus >( "1", 20, true ) ) );

    BOOST_CHECK( not breakpointLocusPtrComp()( std::make_shared< BreakpointLocus >( "1", 20, true ),
                                               std::make_shared< BreakpointLocus >( "1", 20, false ) ) );

    BOOST_CHECK( not breakpointLocusPtrComp()( std::make_shared< BreakpointLocus >( "1", 20, false ),
                                               std::make_shared< BreakpointLocus >( "1", 20, false ) ) );

    BOOST_CHECK( not breakpointLocusPtrComp()( std::make_shared< BreakpointLocus >( "1", 20, true ),
                                               std::make_shared< BreakpointLocus >( "1", 20, true ) ) );
}

BOOST_AUTO_TEST_CASE( testStorageOfMateRegions )
{
    BreakpointLocus breakpoint( "1", 1, true );
    breakpoint.addMateRegion( Region( "1", 0, 100 ), 0 );
    breakpoint.addMateRegion( Region( "2", 0, 100 ), 0 );

    const auto mateRegions = breakpoint.mateRegions();
    std::set< Region > expected = {Region( "1", 0, 100 ), Region( "2", 0, 100 )};

    BOOST_CHECK_EQUAL_COLLECTIONS( mateRegions.cbegin(), mateRegions.cend(), expected.cbegin(), expected.cend() );
}

BOOST_AUTO_TEST_CASE( testShouldMergeOverlappingRegions )
{

    BreakpointLocus breakpoint( "1", 1, true );
    breakpoint.addMateRegion( Region( "1", 0, 100 ), 0 );
    breakpoint.addMateRegion( Region( "1", 99, 200 ), 0 );

    const auto mateRegions = breakpoint.mateRegions();
    std::set< Region > expected = {Region( "1", 0, 200 )};

    BOOST_CHECK_EQUAL_COLLECTIONS( mateRegions.cbegin(), mateRegions.cend(), expected.cbegin(), expected.cend() );
}

BOOST_AUTO_TEST_CASE( testShouldMergeThreeOverlappingRegions )
{

    BreakpointLocus breakpoint( "1", 1, true );
    breakpoint.addMateRegion( Region( "1", 0, 100 ), 0 );
    breakpoint.addMateRegion( Region( "1", 150, 200 ), 0 );
    breakpoint.addMateRegion( Region( "1", 50, 160 ), 0 );

    const auto mateRegions = breakpoint.mateRegions();
    std::set< Region > expected = {Region( "1", 0, 200 )};

    BOOST_CHECK_EQUAL_COLLECTIONS( mateRegions.cbegin(), mateRegions.cend(), expected.cbegin(), expected.cend() );
}
