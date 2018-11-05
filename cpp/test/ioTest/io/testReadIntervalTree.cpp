// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "io/read.hpp"
#include "io/readIntervalTree.hpp"
#include "io/readRange.hpp"

using Read = wecall::io::Read;
using Cigar = wecall::alignment::Cigar;
using Interval = wecall::utils::Interval;
using wecall::utils::BasePairSequence;
using wecall::io::RegionsReads;
using wecall::caller::SetRegions;
using wecall::caller::Region;

BOOST_AUTO_TEST_CASE( shouldFilterReadsWithLowMappingQuality )
{
    wecall::io::readIntervalTree_t readContainer( 0, 100 );

    std::string qname = "test";
    auto refSequence =
        std::make_shared< wecall::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );

    uint8_t mappingQuality = 10;
    readContainer.insert( std::make_shared< Read >( std::string( 2, 'A' ), std::string( 2, 'Q' ), qname, Cigar( "2M" ),
                                                    0, 1, 0, mappingQuality, 0, 0, 0, refSequence ) );

    SetRegions setRegions( refSequence->region() );
    const auto readRange = readContainer.getFullRange();

    {
        const RegionsReads range( setRegions, readRange, mappingQuality );
        BOOST_CHECK_EQUAL( std::distance( range.begin(), range.end() ), 1 );

        const auto subReads = range.getSubRegionReads( setRegions );
        BOOST_CHECK_EQUAL( std::distance( subReads.begin(), subReads.end() ), 1 );
    }
    {
        const RegionsReads range( setRegions, readRange, mappingQuality + 1 );
        BOOST_CHECK_EQUAL( std::distance( range.begin(), range.end() ), 0 );

        const auto subReads = range.getSubRegionReads( setRegions );
        BOOST_CHECK_EQUAL( std::distance( subReads.begin(), subReads.end() ), 0 );
    }
}

BOOST_AUTO_TEST_CASE( shouldGetCorrectSubrangesOneReadWithMatches )
{
    wecall::io::readIntervalTree_t readContainer( 0, 100 );
    std::string qname = "test";
    auto refSequence =
        std::make_shared< wecall::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );

    readContainer.insert( std::make_shared< Read >( std::string( 2, 'A' ), std::string( 2, 'Q' ), qname, Cigar( "2M" ),
                                                    0, 1, 0, 0, 0, 0, 0, refSequence ) );

    auto range0 = readContainer.getFullRange();
    BOOST_CHECK_EQUAL( std::distance( range0.first, range0.second ), 1 );

    auto range1 = readContainer.getSubRange( Interval( 0, 1 ) );
    BOOST_CHECK_EQUAL( std::distance( range1.first, range1.second ), 0 );

    auto range2 = readContainer.getSubRange( Interval( 1, 2 ) );
    BOOST_CHECK_EQUAL( std::distance( range2.first, range2.second ), 1 );

    auto range3 = readContainer.getSubRange( Interval( 2, 9 ) );
    BOOST_CHECK_EQUAL( std::distance( range3.first, range3.second ), 1 );

    auto range4 = readContainer.getSubRange( Interval( 3, 10 ) );
    BOOST_CHECK_EQUAL( std::distance( range4.first, range4.second ), 0 );

    for ( const auto & read : readContainer )
    {
        BOOST_CHECK_EQUAL( read.getReadGroupID(), qname );
    }
}

BOOST_AUTO_TEST_CASE( shouldGetCorrectSubrangesOneReadWithPureInsertion )
{
    auto refSequence =
        std::make_shared< wecall::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );
    wecall::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( std::make_shared< Read >( std::string( 0, 'A' ), std::string( 0, 'Q' ), "0", Cigar( "0I" ), 0,
                                                    1, 0, 0, 0, 0, 0, refSequence ) );

    auto range0 = readContainer.getFullRange();
    BOOST_CHECK_EQUAL( std::distance( range0.first, range0.second ), 1 );

    auto range1 = readContainer.getSubRange( Interval( 0, 1 ) );
    BOOST_CHECK_EQUAL( std::distance( range1.first, range1.second ), 0 );

    auto range2 = readContainer.getSubRange( Interval( 1, 2 ) );
    BOOST_CHECK_EQUAL( std::distance( range2.first, range2.second ), 0 );

    auto range3 = readContainer.getSubRange( Interval( 0, 2 ) );
    BOOST_CHECK_EQUAL( std::distance( range3.first, range3.second ), 1 );

    auto range4 = readContainer.getSubRange( Interval( 1, 1 ) );
    BOOST_CHECK_EQUAL( std::distance( range4.first, range4.second ), 1 );

    auto range5 = readContainer.getSubRange( Interval( 2, 3 ) );
    BOOST_CHECK_EQUAL( std::distance( range5.first, range5.second ), 0 );
}

BOOST_AUTO_TEST_CASE( shouldGetCorrectSubrangesOneReadWithInsertionAtStart )
{
    auto refSequence =
        std::make_shared< wecall::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );
    wecall::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( std::make_shared< Read >( std::string( 8, 'A' ), std::string( 8, 'Q' ), "0", Cigar( "7I1M" ),
                                                    0, 1, 0, 0, 0, 0, 0, refSequence ) );

    auto range0 = readContainer.getFullRange();
    BOOST_CHECK_EQUAL( std::distance( range0.first, range0.second ), 1 );

    auto range1 = readContainer.getSubRange( Interval( 0, 1 ) );
    BOOST_CHECK_EQUAL( std::distance( range1.first, range1.second ), 1 );

    auto range2 = readContainer.getSubRange( Interval( 1, 2 ) );
    BOOST_CHECK_EQUAL( std::distance( range2.first, range2.second ), 1 );

    auto range3 = readContainer.getSubRange( Interval( 0, 2 ) );
    BOOST_CHECK_EQUAL( std::distance( range3.first, range3.second ), 1 );

    auto range4 = readContainer.getSubRange( Interval( 1, 1 ) );
    BOOST_CHECK_EQUAL( std::distance( range4.first, range4.second ), 1 );

    auto range5 = readContainer.getSubRange( Interval( 2, 3 ) );
    BOOST_CHECK_EQUAL( std::distance( range5.first, range5.second ), 0 );
}

BOOST_AUTO_TEST_CASE( testRegionSetReadsIteration )
{
    auto refSequence =
        std::make_shared< wecall::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );
    wecall::io::readIntervalTree_t readContainer( 0, 100 );
    const int64_t startPos = 1;
    readContainer.insert( std::make_shared< Read >( BasePairSequence( 8, 'A' ), std::string( 8, 'Q' ), "0",
                                                    Cigar( "8M" ), 0, startPos, 0, 0, 0, 0, 0, refSequence ) );

    const auto readRange = readContainer.getFullRange();
    SetRegions setRegions;
    setRegions.insert( Region( "1", 1, 2 ) );

    RegionsReads regionSetReads( setRegions, readRange, 0 );
    for ( const auto & read : regionSetReads )
    {
        BOOST_CHECK_EQUAL( read.sequence(), BasePairSequence( 8, 'A' ) );
    }
    BOOST_CHECK_EQUAL( 1, std::distance( regionSetReads.begin(), regionSetReads.end() ) );
}

BOOST_AUTO_TEST_CASE( testShouldOnlyRetrieveReadsThatOverlapARegion )
{
    wecall::io::readIntervalTree_t readContainer( 0, 100 );
    const int64_t startPos1 = 1;
    const std::size_t length = 8;
    auto refSequence =
        std::make_shared< wecall::utils::ReferenceSequence >( Region( "1", 0, 20 ), std::string( 20, 'A' ) );
    readContainer.insert( std::make_shared< Read >( BasePairSequence( length, 'A' ), std::string( length, 'Q' ), "0",
                                                    Cigar( std::to_string( length ) + "M" ), 0, startPos1, 0, 0, 0, 0,
                                                    0, refSequence ) );
    const int64_t startPos2 = 10;
    readContainer.insert( std::make_shared< Read >( BasePairSequence( length, 'A' ), std::string( length, 'Q' ), "0",
                                                    Cigar( std::to_string( length ) + "M" ), 0, startPos2, 0, 0, 0, 0,
                                                    0, refSequence ) );

    const auto readRange = readContainer.getFullRange();
    SetRegions setRegions;
    setRegions.insert( Region( "1", 0, 1 ) );
    setRegions.insert( Region( "1", 9, 10 ) );
    setRegions.insert( Region( "1", 18, 20 ) );

    RegionsReads regionSetReads( setRegions, readRange, 0 );

    BOOST_CHECK_EQUAL( 0, std::distance( regionSetReads.begin(), regionSetReads.end() ) );
}

BOOST_AUTO_TEST_CASE( testGetSubRegionReadsShouldThrowIfSubRegionIsNotInContained )
{
    wecall::io::readIntervalTree_t readContainer( 0, 100 );
    const int64_t startPos1 = 1;
    const std::size_t length = 8;
    auto refSequence =
        std::make_shared< wecall::utils::ReferenceSequence >( Region( "1", 0, 20 ), std::string( 20, 'A' ) );
    readContainer.insert( std::make_shared< Read >( BasePairSequence( length, 'A' ), std::string( length, 'Q' ), "0",
                                                    Cigar( std::to_string( length ) + "M" ), 0, startPos1, 0, 0, 0, 0,
                                                    0, refSequence ) );
    const int64_t startPos2 = 10;
    readContainer.insert( std::make_shared< Read >( BasePairSequence( length, 'A' ), std::string( length, 'Q' ), "0",
                                                    Cigar( std::to_string( length ) + "M" ), 0, startPos2, 0, 0, 0, 0,
                                                    0, refSequence ) );

    const auto readRange = readContainer.getFullRange();
    SetRegions setRegions;
    setRegions.insert( Region( "1", 9, 10 ) );
    setRegions.insert( Region( "1", 18, 20 ) );

    SetRegions setRegions2;
    setRegions2.insert( Region( "1", 8, 10 ) );
    setRegions2.insert( Region( "1", 18, 20 ) );

    RegionsReads regionSetReads( setRegions, readRange, 0 );

    BOOST_CHECK_THROW( regionSetReads.getSubRegionReads( setRegions2 ), wecall::utils::wecall_exception );
}
