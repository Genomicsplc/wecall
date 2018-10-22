// All content Copyright (C) 2018 Genomics plc
#include "io/readDataSet.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

using wecall::alignment::Cigar;
using wecall::caller::Region;
using wecall::io::Read;
using wecall::io::ReadDataset;
using wecall::io::RegionsReads;
using wecall::utils::BasePairSequence;

BOOST_AUTO_TEST_CASE( notOverlappingForReadOnlyOverlappingRegion1 )
{
    Region region = Region( "1", 0, 30 );
    Region region1 = Region( "1", 10, 11 );
    Region region2 = Region( "1", 20, 21 );
    auto refSequence = std::make_shared< wecall::utils::ReferenceSequence >( region, std::string( 30, 'A' ) );

    const int64_t startPos = 10;
    const uint8_t mapQual = 10;
    const auto read1 = std::make_shared< Read >( BasePairSequence( std::string( 10, 'C' ) ), std::string( 10, 'Q' ),
                                                 "0", Cigar( "10M" ), 0, startPos, 0, mapQual, 0, 0, 0, refSequence );

    wecall::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );

    RegionsReads regionSetReads( region, readContainer.getFullRange(), 0 );

    wecall::io::perSampleRegionsReads_t reads = {{"sample1", regionSetReads}};

    bool overlaps = wecall::io::readsOverlappingRegions( reads, region1, region2 );

    BOOST_CHECK( not overlaps );
}

BOOST_AUTO_TEST_CASE( notOverlappingForReadOnlyOverlappingRegion2 )
{
    std::cout << "running test" << std::endl;
    Region region = Region( "1", 0, 30 );
    Region region1 = Region( "1", 10, 11 );
    Region region2 = Region( "1", 20, 21 );
    auto refSequence = std::make_shared< wecall::utils::ReferenceSequence >( region, std::string( 30, 'A' ) );

    const int64_t startPos = 11;
    const uint8_t mapQual = 10;
    const auto read1 = std::make_shared< Read >( BasePairSequence( std::string( 10, 'C' ) ), std::string( 10, 'Q' ),
                                                 "0", Cigar( "10M" ), 0, startPos, 0, mapQual, 0, 0, 0, refSequence );

    wecall::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );

    RegionsReads regionSetReads( region, readContainer.getFullRange(), 0 );

    wecall::io::perSampleRegionsReads_t reads = {{"sample1", regionSetReads}};

    bool overlaps = wecall::io::readsOverlappingRegions( reads, region1, region2 );

    BOOST_CHECK( not overlaps );
}

BOOST_AUTO_TEST_CASE( overlappingForReadJustOverlapping )
{
    Region region = Region( "1", 0, 30 );
    Region region1 = Region( "1", 10, 11 );
    Region region2 = Region( "1", 20, 21 );
    auto refSequence = std::make_shared< wecall::utils::ReferenceSequence >( region, std::string( 30, 'A' ) );

    const int64_t startPos = 10;
    const uint8_t mapQual = 10;
    const auto read1 = std::make_shared< Read >( BasePairSequence( std::string( 11, 'C' ) ), std::string( 11, 'Q' ),
                                                 "0", Cigar( "11M" ), 0, startPos, 0, mapQual, 0, 0, 0, refSequence );

    wecall::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );

    RegionsReads regionSetReads( region, readContainer.getFullRange(), 0 );

    wecall::io::perSampleRegionsReads_t reads = {{"sample1", regionSetReads}};

    bool overlaps = wecall::io::readsOverlappingRegions( reads, region1, region2 );

    BOOST_CHECK( overlaps );
}

BOOST_AUTO_TEST_CASE( notOverlappingForReadOutsideOfRange )
{
    Region region = Region( "1", 0, 30 );
    Region region1 = Region( "1", 10, 11 );
    Region region2 = Region( "1", 20, 21 );
    auto refSequence = std::make_shared< wecall::utils::ReferenceSequence >( region, std::string( 30, 'A' ) );

    const int64_t startPos = 0;
    const uint8_t mapQual = 10;
    const auto read1 = std::make_shared< Read >( BasePairSequence( std::string( 5, 'C' ) ), std::string( 5, 'Q' ), "0",
                                                 Cigar( "5M" ), 0, startPos, 0, mapQual, 0, 0, 0, refSequence );

    wecall::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );

    RegionsReads regionSetReads( region, readContainer.getFullRange(), 0 );

    wecall::io::perSampleRegionsReads_t reads = {{"sample1", regionSetReads}};

    bool overlaps = wecall::io::readsOverlappingRegions( reads, region1, region2 );

    BOOST_CHECK( not overlaps );
}

BOOST_AUTO_TEST_CASE( overlappingForTwoSamples )
{
    Region region = Region( "1", 0, 30 );
    Region region1 = Region( "1", 10, 11 );
    Region region2 = Region( "1", 20, 21 );
    auto refSequence = std::make_shared< wecall::utils::ReferenceSequence >( region, std::string( 30, 'A' ) );

    const int64_t startPos = 10;
    const uint8_t mapQual = 10;
    const auto overlappingRead =
        std::make_shared< Read >( BasePairSequence( std::string( 11, 'C' ) ), std::string( 11, 'Q' ), "0",
                                  Cigar( "11M" ), 0, startPos, 0, mapQual, 0, 0, 0, refSequence );

    const auto nonOverlappingRead =
        std::make_shared< Read >( BasePairSequence( std::string( 5, 'C' ) ), std::string( 5, 'Q' ), "0", Cigar( "5M" ),
                                  0, startPos, 0, mapQual, 0, 0, 0, refSequence );

    wecall::io::readIntervalTree_t readContainer1( 0, 100 );
    readContainer1.insert( nonOverlappingRead );

    wecall::io::readIntervalTree_t readContainer2( 0, 100 );
    readContainer2.insert( overlappingRead );

    RegionsReads regionSetReads1( region, readContainer1.getFullRange(), 0 );
    RegionsReads regionSetReads2( region, readContainer2.getFullRange(), 0 );

    wecall::io::perSampleRegionsReads_t reads = {{"sample1", regionSetReads1}, {"sample2", regionSetReads2}};

    bool overlaps = wecall::io::readsOverlappingRegions( reads, region1, region2 );

    BOOST_CHECK( overlaps );
}

BOOST_AUTO_TEST_CASE( nonOverlappingForTwoSamples )
{
    Region region = Region( "1", 0, 30 );
    Region region1 = Region( "1", 10, 11 );
    Region region2 = Region( "1", 20, 21 );
    auto refSequence = std::make_shared< wecall::utils::ReferenceSequence >( region, std::string( 30, 'A' ) );

    const int64_t startPos = 10;
    const uint8_t mapQual = 10;
    const auto nonOverlappingRead1 =
        std::make_shared< Read >( BasePairSequence( std::string( 10, 'C' ) ), std::string( 10, 'Q' ), "0",
                                  Cigar( "10M" ), 0, startPos, 0, mapQual, 0, 0, 0, refSequence );

    const auto nonOverlappingRead2 =
        std::make_shared< Read >( BasePairSequence( std::string( 5, 'C' ) ), std::string( 5, 'Q' ), "0", Cigar( "5M" ),
                                  0, startPos, 0, mapQual, 0, 0, 0, refSequence );

    wecall::io::readIntervalTree_t readContainer1( 0, 100 );
    readContainer1.insert( nonOverlappingRead1 );

    wecall::io::readIntervalTree_t readContainer2( 0, 100 );
    readContainer2.insert( nonOverlappingRead2 );

    RegionsReads regionSetReads1( region, readContainer1.getFullRange(), 0 );
    RegionsReads regionSetReads2( region, readContainer2.getFullRange(), 0 );

    wecall::io::perSampleRegionsReads_t reads = {{"sample1", regionSetReads1}, {"sample2", regionSetReads2}};

    bool overlaps = wecall::io::readsOverlappingRegions( reads, region1, region2 );

    BOOST_CHECK( not overlaps );
}