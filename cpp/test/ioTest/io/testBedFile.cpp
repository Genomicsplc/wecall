#include "io/bamFile.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "ioFixture.hpp"

using echidna::caller::Region;
using echidna::caller::regions_t;

using namespace echidna::test;

BOOST_FIXTURE_TEST_CASE( shouldReadBedWithSingleThreeColumnLine, BedFileFixture )
{
    regions_t expectedRegions = {Region( "chr1", 1, 2 )};
    auto regions = bedFiles.at( 0 )->getRegions();
    BOOST_CHECK_EQUAL( regions.size(), expectedRegions.size() );
    BOOST_CHECK_EQUAL_COLLECTIONS( regions.begin(), regions.end(), expectedRegions.begin(), expectedRegions.end() );
}

BOOST_FIXTURE_TEST_CASE( shouldReadBedWithMultipleThreeColumnLines, BedFileFixture )
{
    regions_t expectedRegions = {Region( "chr1", 2, 3 ), Region( "chr2", 5, 7 )};
    auto regions = bedFiles.at( 1 )->getRegions();
    BOOST_CHECK_EQUAL( regions.size(), expectedRegions.size() );
    BOOST_CHECK_EQUAL_COLLECTIONS( regions.begin(), regions.end(), expectedRegions.begin(), expectedRegions.end() );
}

BOOST_FIXTURE_TEST_CASE( shouldReadRealExampleWithoutHeader, BedFileFixture )
{
    // example take from https://genome.ucsc.edu/FAQ/FAQformat.html#format1
    regions_t expectedRegions = {Region( "chr7", 0, 10 ), Region( "chr7", 20, 25 ), Region( "chr7", 30, 50 )};
    auto regions = bedFiles.at( 2 )->getRegions();
    BOOST_CHECK_EQUAL( regions.size(), expectedRegions.size() );
    BOOST_CHECK_EQUAL_COLLECTIONS( regions.begin(), regions.end(), expectedRegions.begin(), expectedRegions.end() );
}

BOOST_FIXTURE_TEST_CASE( shouldReadRealExampleWithHeader, BedFileFixture )
{
    // example take from https://genome.ucsc.edu/FAQ/FAQformat.html#format1
    regions_t expectedRegions = {
        Region( "chr7", 0, 10 ), Region( "chr7", 20, 25 ), Region( "gi|734691289|gb|KN707645.1|", 30, 50 )};
    auto regions = bedFiles.at( 3 )->getRegions();
    BOOST_CHECK_EQUAL( regions.size(), expectedRegions.size() );
    BOOST_CHECK_EQUAL_COLLECTIONS( regions.begin(), regions.end(), expectedRegions.begin(), expectedRegions.end() );
}
