// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "caller/region.hpp"
#include "caller/regionUtils.hpp"
#include "utils/flatten.hpp"

using wecall::caller::Region;
using wecall::utils::Interval;
using wecall::caller::parseRegionString;

BOOST_AUTO_TEST_CASE( testGettingContigsFromRegions )
{
    std::vector< Region > regions = {{"1", 1, 2}, {"1", 5, 9}, {"2", 0, 20}};
    const std::map< std::string, Interval > referenceContigs = {{"1", Interval( 0, 100 )}, {"2", Interval( 0, 100 )}};

    std::vector< Region > expectedContigs = {Region( "1", Interval( 0, 100 ) ), Region( "2", Interval( 0, 100 ) )};

    auto contigs = wecall::caller::getContigsFromRegions( regions, referenceContigs );

    BOOST_CHECK_EQUAL( contigs.size(), 2 );
    BOOST_CHECK_EQUAL( contigs.front().contig(), expectedContigs.front().contig() );  // Doesn't change order
    BOOST_CHECK_EQUAL( contigs.back().contig(), expectedContigs.back().contig() );    // Doesn't change order
    BOOST_CHECK_EQUAL( contigs.front().start(), expectedContigs.front().start() );    // Doesn't change order
    BOOST_CHECK_EQUAL( contigs.back().start(), expectedContigs.back().start() );      // Doesn't change order
    BOOST_CHECK_EQUAL( contigs.front().end(), expectedContigs.front().end() );        // Doesn't change order
    BOOST_CHECK_EQUAL( contigs.back().end(), expectedContigs.back().end() );          // Doesn't change order
}

BOOST_AUTO_TEST_CASE( padRegionsShouldReturnEmptyVectorForEmptyInput )
{
    std::vector< Region > input = {};
    std::vector< Region > expectedOutput = {};

    auto results = wecall::caller::padRegions( input, 0 );
    BOOST_CHECK_EQUAL_COLLECTIONS( results.begin(), results.end(), expectedOutput.begin(), expectedOutput.end() );
}

BOOST_AUTO_TEST_CASE( mergeRegionsShouldReturnEmptyVectorForEmptyInput )
{
    std::vector< Region > input = {};
    std::vector< Region > expectedOutput = {};

    auto results = wecall::caller::mergeRegions( input );
    BOOST_CHECK_EQUAL_COLLECTIONS( results.begin(), results.end(), expectedOutput.begin(), expectedOutput.end() );
}

BOOST_AUTO_TEST_CASE( padRegionsShouldFailIfNotAllContigsAreTheSame )
{
    std::vector< Region > input = {{"1", 1, 2}, {"2", 1, 2}};

    BOOST_CHECK_THROW( wecall::caller::padRegions( input, 0 ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( mergeRegionsShouldFailIfNotAllContigsAreTheSame )
{
    std::vector< Region > input = {{"1", 1, 2}, {"2", 1, 2}};

    BOOST_CHECK_THROW( wecall::caller::mergeRegions( input ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( padRegionsShouldFailIfNotAllRegionsAreSorted )
{
    std::vector< Region > input = {{"1", 1, 2}, {"1", 0, 20}};

    BOOST_CHECK_THROW( wecall::caller::padRegions( input, 0 ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( mergeRegionsShouldFailIfNotAllRegionsAreSorted )
{
    std::vector< Region > input = {{"1", 1, 2}, {"1", 0, 20}};

    BOOST_CHECK_THROW( wecall::caller::mergeRegions( input ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( padRegionsShouldPadRegionsIntoSingleRegion )
{
    std::vector< Region > input = {{"1", 1, 2}, {"1", 3, 4}};
    const auto padding = 1;
    std::vector< Region > expectedOutput = {{"1", 0, 5}};

    auto results = wecall::caller::padRegions( input, padding );
    BOOST_CHECK_EQUAL_COLLECTIONS( results.begin(), results.end(), expectedOutput.begin(), expectedOutput.end() );
}

BOOST_AUTO_TEST_CASE( padParitionedShouldMaintainContigs )
{
    std::vector< std::vector< Region > > input = {{{"1", 1, 2}, {"1", 3, 4}}, {{"2", 1, 2}}};
    const auto padding = 1;
    std::vector< std::vector< Region > > expectedOutput = {{{"1", 0, 5}}, {{"2", 0, 3}}};

    auto results = wecall::caller::padPartitionedRegions( input, padding );
    BOOST_REQUIRE_EQUAL( results.size(), 2 );

    BOOST_CHECK_EQUAL_COLLECTIONS( results[0].begin(), results[0].end(), expectedOutput[0].begin(),
                                   expectedOutput[0].end() );
    BOOST_CHECK_EQUAL_COLLECTIONS( results[1].begin(), results[1].end(), expectedOutput[1].begin(),
                                   expectedOutput[1].end() );
}

BOOST_AUTO_TEST_CASE( mergeRegionsShouldNotMergeNonOverlappingRegions )
{
    std::vector< Region > input = {{"1", 1, 2}, {"1", 3, 4}};
    std::vector< Region > expectedOutput = {{"1", 1, 2}, {"1", 3, 4}};

    auto results = wecall::caller::mergeRegions( input );
    BOOST_CHECK_EQUAL_COLLECTIONS( results.begin(), results.end(), expectedOutput.begin(), expectedOutput.end() );
}

BOOST_AUTO_TEST_CASE( mergeRegionsShouldMergeOverlappingRegions )
{
    std::vector< Region > input = {{"1", 1, 4}, {"1", 3, 5}};
    std::vector< Region > expectedOutput = {{"1", 1, 5}};

    auto results = wecall::caller::mergeRegions( input );
    BOOST_CHECK_EQUAL_COLLECTIONS( results.begin(), results.end(), expectedOutput.begin(), expectedOutput.end() );
}

BOOST_AUTO_TEST_CASE( mergeRegionsShouldMergeTouchingRegions )
{
    std::vector< Region > input = {{"1", 1, 3}, {"1", 3, 5}};
    std::vector< Region > expectedOutput = {{"1", 1, 5}};

    auto results = wecall::caller::mergeRegions( input );
    BOOST_CHECK_EQUAL_COLLECTIONS( results.begin(), results.end(), expectedOutput.begin(), expectedOutput.end() );
}

BOOST_AUTO_TEST_CASE( testParseRegionStringJustContig )
{
    std::string contig = "1";
    const std::map< std::string, Interval > referenceContigs = {{"1", Interval( 0, 0 )}};

    Region region = parseRegionString( contig, referenceContigs );

    BOOST_CHECK_EQUAL( region.contig(), contig );
    BOOST_CHECK_EQUAL( region.start(), 0 );
    BOOST_CHECK_EQUAL( region.end(), 0 );
    BOOST_CHECK( region.hasNoRange() );
    BOOST_CHECK_EQUAL( region.toString(), "1" );
}
BOOST_AUTO_TEST_CASE( testParseRegionStringJustContigWithReference )
{
    std::string contig = "1";
    const std::map< std::string, Interval > referenceContigs = {{"1", Interval( 10, 100 )}};

    Region region = parseRegionString( contig, referenceContigs );

    BOOST_CHECK_EQUAL( region.contig(), contig );
    BOOST_CHECK_EQUAL( region.start(), 10 );
    BOOST_CHECK_EQUAL( region.end(), 100 );
    BOOST_CHECK( not region.hasNoRange() );
    BOOST_CHECK_EQUAL( region.toString(), "1:10-100" );
}

BOOST_AUTO_TEST_CASE( testParseFullWithReference )
{
    std::string strRegion = "1:0-10";
    const std::map< std::string, Interval > referenceContigs = {{"1", Interval( 0, 100 )}};

    Region region = parseRegionString( strRegion, referenceContigs );

    BOOST_CHECK_EQUAL( region.contig(), "1" );
    BOOST_CHECK_EQUAL( region.start(), 0 );
    BOOST_CHECK_EQUAL( region.end(), 10 );
    BOOST_CHECK( region.end() != referenceContigs.at( "1" ).end() );
    BOOST_CHECK( not region.hasNoRange() );
    BOOST_CHECK_EQUAL( region.toString(), strRegion );
}

BOOST_AUTO_TEST_CASE( testParseRegionStringComplexContig )
{
    std::string contig = "gi|725026149|gb|JTFH01000117.1|";
    const std::map< std::string, Interval > referenceContigs = {{contig, Interval( 0, 0 )}};

    Region region = parseRegionString( contig, referenceContigs );

    BOOST_CHECK_EQUAL( region.contig(), contig );
    BOOST_CHECK_EQUAL( region.start(), 0 );
    BOOST_CHECK_EQUAL( region.end(), 0 );
    BOOST_CHECK( region.hasNoRange() );
    BOOST_CHECK_EQUAL( region.toString(), contig );
}

BOOST_AUTO_TEST_CASE( testParseRegionContigStartPos )
{
    const std::string contig = "gi|725026149|gb|JTFH01000117.1|";

    std::string regionString = contig + ":1-13";
    const std::map< std::string, Interval > referenceContigs = {{contig, Interval( 0, 100 )}};
    Region region = parseRegionString( regionString, referenceContigs );

    BOOST_CHECK_EQUAL( region.contig(), "gi|725026149|gb|JTFH01000117.1|" );
    BOOST_CHECK_EQUAL( region.start(), 1 );
    BOOST_CHECK_EQUAL( region.end(), 13 );
    BOOST_CHECK( not region.hasNoRange() );
    BOOST_CHECK_EQUAL( region.toString(), regionString );
    BOOST_CHECK_EQUAL( region.size(), 13 - 1 );
    BOOST_CHECK_EQUAL( region.toString(), regionString );
}

BOOST_AUTO_TEST_CASE( testFailsIfNoDashWithColon )
{
    const std::map< std::string, Interval > contigs = {{"2", Interval( 0, 100 )}};
    BOOST_CHECK_THROW( parseRegionString( "2:123", contigs ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( testFailsIfTwoColons )
{
    const std::map< std::string, Interval > contigs = {{"1", Interval( 0, 100 )}, {"2", Interval( 0, 100 )}};
    BOOST_CHECK_THROW( parseRegionString( "2:1:23-26", contigs ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( testFailsIfHaveNonNumericRanges )
{
    const std::map< std::string, Interval > contigs = {{"1", Interval( 0, 100 )}, {"2", Interval( 0, 100 )}};
    BOOST_CHECK_THROW( parseRegionString( "2:Inf-12", contigs ), wecall::utils::wecall_exception );
    BOOST_CHECK_THROW( parseRegionString( "2:1-Inf", contigs ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( testFailsIfTooManyDashes )
{
    const std::map< std::string, Interval > contigs = {{"1", Interval( 0, 100 )}, {"2", Interval( 0, 100 )}};
    BOOST_CHECK_THROW( parseRegionString( "2:A-12-14", contigs ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( testFailsIfNotValidHalfOpenRegion )
{
    const std::map< std::string, Interval > contigs = {{"1", Interval( 0, 100 )}, {"2", Interval( 0, 100 )}};
    BOOST_CHECK_THROW( parseRegionString( "2:10-9", contigs ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( testFailsIfContainsNegativeValues )
{
    const std::map< std::string, Interval > contigs = {{"1", Interval( 0, 100 )}, {"2", Interval( 0, 100 )}};
    BOOST_CHECK_THROW( parseRegionString( "2:-10-10", contigs ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( testSetsTheIntervalToBeEmptyAtZeroIfContigNotContainedInReferenceContigs )
{
    const std::map< std::string, Interval > contigs = {{"1", Interval( 0, 100 )}, {"2", Interval( 0, 100 )}};
    BOOST_CHECK_EQUAL( parseRegionString( "3", contigs ), Region( "3", 0, 0 ) );
}
