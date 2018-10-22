// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "readrecalibration/readRecalibration.hpp"
#include "io/read.hpp"
#include "alignment/cigar.hpp"

#include <memory>

using Read = wecall::io::Read;
using wecall::io::RegionsReads;
using wecall::caller::Region;
using wecall::alignment::Cigar;
using wecall::alignment::cigarFlags;

BOOST_AUTO_TEST_CASE( testFloorLowQualityScores )
{
    char lowQual = 5;

    auto refSequence =
        std::make_shared< wecall::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );

    const auto read1 = std::make_shared< Read >( std::string( 3, 'T' ), std::string( 3, lowQual ), "", Cigar( "3M" ), 0,
                                                 0, 0, 1, 2, 2, 0, refSequence, "read1" );

    wecall::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );

    RegionsReads regionSetReads( Region( "1", 1, 2 ), readContainer.getFullRange(), 0 );
    wecall::io::perSampleRegionsReads_t perSampleReads = {{"sample1", regionSetReads}};

    wecall::corrector::floorLowQualityScores( perSampleReads, 5, 2 );

    const auto qualities = readContainer.begin()->getQualities();
    BOOST_REQUIRE_EQUAL( qualities.size(), 3 );
    for ( const auto & qualityChar : qualities )
    {
        BOOST_CHECK_EQUAL( qualityChar, 2 );
    }
}

BOOST_AUTO_TEST_CASE( testShouldNotFloorHighQualityScores )
{
    char highQual = 5;

    auto refSequence =
        std::make_shared< wecall::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );

    const auto read1 = std::make_shared< Read >( std::string( 3, 'T' ), std::string( 3, highQual ), "", Cigar( "3M" ),
                                                 0, 0, 0, 1, 2, 2, 0, refSequence, "read1" );

    wecall::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );

    RegionsReads regionSetReads( Region( "1", 1, 2 ), readContainer.getFullRange(), 0 );
    wecall::io::perSampleRegionsReads_t perSampleReads = {{"sample1", regionSetReads}};

    const auto qualities = readContainer.begin()->getQualities();
    BOOST_REQUIRE_EQUAL( qualities.size(), 3 );
    for ( const auto & qualityChar : qualities )
    {
        BOOST_CHECK_EQUAL( qualityChar, highQual );
    }
}
