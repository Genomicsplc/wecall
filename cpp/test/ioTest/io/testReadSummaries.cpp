// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "io/read.hpp"
#include "io/readIntervalTree.hpp"
#include "io/readRange.hpp"
#include "io/readSummaries.hpp"

using Read = wecall::io::Read;
using Cigar = wecall::alignment::Cigar;
using wecall::caller::SetRegions;
using wecall::caller::Region;
using wecall::io::RegionsReads;
using wecall::utils::BasePairSequence;
using readCoverage_t = wecall::io::readsummaries::readCoverage_t;

BOOST_AUTO_TEST_CASE( shouldComputeCoverageDelatsFor1Read )
{
    auto refSequence =
        std::make_shared< wecall::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'A' ) );

    const int64_t startPos = 0;
    const auto read1 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    wecall::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );

    RegionsReads regionSetReads( Region( "1", 1, 2 ), readContainer.getFullRange(), 0 );

    wecall::io::perSampleRegionsReads_t perSampleReads = {{"sample1", regionSetReads}};
    Region subRegion = Region( "1", 1, 3 );

    const auto readCoverageDeltas =
        wecall::io::readsummaries::getReadCoverageDeltas( perSampleReads, subRegion, 1, subRegion.size() );

    BOOST_REQUIRE_EQUAL( readCoverageDeltas.size(), 1 );
    BOOST_REQUIRE_EQUAL( readCoverageDeltas[0].size(), 2 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][0], 1 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][1], 0 );  // no change to base level
}

BOOST_AUTO_TEST_CASE( shouldComputeReadDeltasFor1SampleWith2ReadsSpanningWholeRegion )
{
    Region region = Region( "1", 0, 5 );
    auto refSequence = std::make_shared< wecall::utils::ReferenceSequence >( region, std::string( 5, 'A' ) );

    const int64_t startPos = 0;
    const auto read1 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const auto read2 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    wecall::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );
    readContainer.insert( read2 );

    RegionsReads regionSetReads( region, readContainer.getFullRange(), 0 );

    wecall::io::perSampleRegionsReads_t perSampleReads = {{"sample1", regionSetReads}};

    // expected coverage matrix: {{1, 1, 1, 1, 0}, {1, 1, 1, 1, 0};

    const auto readCoverageDeltas =
        wecall::io::readsummaries::getReadCoverageDeltas( perSampleReads, region, 1, region.size() );

    BOOST_REQUIRE_EQUAL( readCoverageDeltas.size(), 1 );
    BOOST_REQUIRE_EQUAL( readCoverageDeltas[0].size(), 5 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][0], 2 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][1], 0 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][2], 0 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][3], 0 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][4], -2 );  // reads end before last position
}

BOOST_AUTO_TEST_CASE( shouldComputeReadDeltasFor2SamplesWith1Read )
{
    Region region = Region( "1", 0, 5 );
    auto refSequence = std::make_shared< wecall::utils::ReferenceSequence >( region, std::string( 5, 'A' ) );

    const int64_t startPos = 0;
    const auto read1 = std::make_shared< Read >( BasePairSequence( 3, 'A' ), std::string( 3, 'Q' ), "0", Cigar( "3M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const auto read2 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos + 1, 0, 0, 0, 0, 0, refSequence );

    wecall::io::readIntervalTree_t readContainer1( 0, 100 );
    readContainer1.insert( read1 );

    wecall::io::readIntervalTree_t readContainer2( 0, 100 );
    readContainer2.insert( read2 );

    RegionsReads regionSetReads1( Region( "1", 1, 2 ), readContainer1.getFullRange(), 0 );
    RegionsReads regionSetReads2( Region( "1", 1, 2 ), readContainer2.getFullRange(), 0 );

    wecall::io::perSampleRegionsReads_t perSampleReads = {{"sample1", regionSetReads1}, {"sample2", regionSetReads2}};

    const auto readCoverageDeltas =
        wecall::io::readsummaries::getReadCoverageDeltas( perSampleReads, region, 2, region.size() );

    // expected coverage matrix: {{1, 1, 1, 0, 0}, {0, 1, 1, 1, 1};

    BOOST_REQUIRE_EQUAL( readCoverageDeltas.size(), 2 );
    BOOST_REQUIRE_EQUAL( readCoverageDeltas[0].size(), 5 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][0], 1 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][1], 0 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][2], 0 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][3], -1 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][4], 0 );  // reads end before last position

    BOOST_CHECK_EQUAL( readCoverageDeltas[1][0], 0 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[1][1], 1 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[1][2], 0 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[1][3], 0 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[1][4], 0 );  // reads end before last position
}

BOOST_AUTO_TEST_CASE( shouldComputeReadDeltasFor2SamplesWith2Reads )
{
    Region region = Region( "1", 0, 5 );
    auto refSequence = std::make_shared< wecall::utils::ReferenceSequence >( region, std::string( 5, 'A' ) );

    const int64_t startPos = 0;
    const auto read1 = std::make_shared< Read >( BasePairSequence( 3, 'A' ), std::string( 3, 'Q' ), "0", Cigar( "3M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const auto read2 = std::make_shared< Read >( BasePairSequence( 3, 'A' ), std::string( 3, 'Q' ), "0", Cigar( "3M" ),
                                                 0, startPos + 1, 0, 0, 0, 0, 0, refSequence );

    const auto read3 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const auto read4 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos + 1, 0, 0, 0, 0, 0, refSequence );

    wecall::io::readIntervalTree_t readContainer1( 0, 100 );
    readContainer1.insert( read1 );
    readContainer1.insert( read2 );

    wecall::io::readIntervalTree_t readContainer2( 0, 100 );
    readContainer2.insert( read3 );
    readContainer2.insert( read4 );

    RegionsReads regionSetReads1( Region( "1", 1, 2 ), readContainer1.getFullRange(), 0 );
    RegionsReads regionSetReads2( Region( "1", 1, 2 ), readContainer2.getFullRange(), 0 );

    wecall::io::perSampleRegionsReads_t perSampleReads = {{"sample1", regionSetReads1}, {"sample2", regionSetReads2}};

    const auto readCoverageDeltas =
        wecall::io::readsummaries::getReadCoverageDeltas( perSampleReads, region, 2, region.size() );

    // expected coverage matrix: {{1, 2, 2, 1, 0}, {1, 2, 2, 2, 1};

    BOOST_REQUIRE_EQUAL( readCoverageDeltas.size(), 2 );
    BOOST_REQUIRE_EQUAL( readCoverageDeltas[0].size(), 5 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][0], 1 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][1], 1 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][2], 0 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][3], -1 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[0][4], -1 );

    BOOST_CHECK_EQUAL( readCoverageDeltas[1][0], 1 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[1][1], 1 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[1][2], 0 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[1][3], 0 );
    BOOST_CHECK_EQUAL( readCoverageDeltas[1][4], -1 );
}

BOOST_AUTO_TEST_CASE( shouldComputeReadDeltaFor1SampleWith3Reads )
{
    const auto region = Region( "1", 0, 10 );
    auto refSequence =
        std::make_shared< wecall::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );

    const int64_t startPos = 0;
    const auto read1 = std::make_shared< Read >( BasePairSequence( 3, 'A' ), std::string( 3, 'Q' ), "0", Cigar( "3M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const auto read2 = std::make_shared< Read >( BasePairSequence( 3, 'A' ), std::string( 3, 'Q' ), "0", Cigar( "3M" ),
                                                 0, startPos + 3, 0, 0, 0, 0, 0, refSequence );

    const auto read3 = std::make_shared< Read >( BasePairSequence( 5, 'A' ), std::string( 5, 'Q' ), "0", Cigar( "5M" ),
                                                 0, startPos + 5, 0, 0, 0, 0, 0, refSequence );

    wecall::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );
    readContainer.insert( read2 );
    readContainer.insert( read3 );

    RegionsReads regionSetReads( region, readContainer.getFullRange(), 0 );

    wecall::io::perSampleRegionsReads_t perSampleReads = {{"sample1", regionSetReads}};

    const auto readCoverageDeltas =
        wecall::io::readsummaries::getReadCoverageDeltas( perSampleReads, region, 1, region.size() );

    // expected coverage matrix: { {1, 1, 1, 1, 1, 2, 1, 1, 1, 1}};
    std::vector< int64_t > exp_deltas( {1, 0, 0, 0, 0, 1, -1, 0, 0, 0} );

    BOOST_REQUIRE_EQUAL( readCoverageDeltas.size(), 1 );
    BOOST_REQUIRE_EQUAL( readCoverageDeltas[0].size(), 10 );
    BOOST_CHECK_EQUAL_COLLECTIONS( readCoverageDeltas[0].begin(), readCoverageDeltas[0].end(), exp_deltas.begin(),
                                   exp_deltas.end() );
}

//-------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( shouldComputeChunkedCoverageFor1SampleWithAlmostConstantReadDepth )
{

    Region region = Region( "1", 0, 5 );

    wecall::io::readsummaries::coverageDeltas_t readDeltas( {{10, 0, 1, 0, 0}} );

    std::vector< readCoverage_t > chunkedReferenceCall =
        wecall::io::readsummaries::getChunkedReferenceCalls( region, 1, region.size(), readDeltas, 0.2 );

    BOOST_REQUIRE_EQUAL( chunkedReferenceCall.size(), 1 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[0].region, region );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[0].minTotalCoverage, 10 );
    BOOST_REQUIRE_EQUAL( chunkedReferenceCall[0].minCoveragePerSample.size(), 1 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[0].minCoveragePerSample[0], 10 );
    BOOST_REQUIRE_EQUAL( chunkedReferenceCall[0].averageCoveragePerSample.size(), 1 );
    BOOST_CHECK_CLOSE( chunkedReferenceCall[0].averageCoveragePerSample[0], 53.0 / 5.0, 1 );
}

BOOST_AUTO_TEST_CASE( shouldComputeChunkedCoverageFor1SampleWithStrongChangesInReadDepth )
{

    Region region = Region( "1", 0, 5 );

    wecall::io::readsummaries::coverageDeltas_t readDeltas( {{5, 0, 2, 0, 0}} );

    std::vector< readCoverage_t > chunkedReferenceCall =
        wecall::io::readsummaries::getChunkedReferenceCalls( region, 1, region.size(), readDeltas, 0.2 );

    BOOST_REQUIRE_EQUAL( chunkedReferenceCall.size(), 2 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[0].region, Region( "1", 0, 2 ) );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[0].minTotalCoverage, 5 );
    BOOST_REQUIRE_EQUAL( chunkedReferenceCall[0].minCoveragePerSample.size(), 1 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[0].minCoveragePerSample[0], 5 );
    BOOST_REQUIRE_EQUAL( chunkedReferenceCall[0].averageCoveragePerSample.size(), 1 );
    BOOST_CHECK_CLOSE( chunkedReferenceCall[0].averageCoveragePerSample[0], 10.0 / 2.0, 1 );

    BOOST_CHECK_EQUAL( chunkedReferenceCall[1].region, Region( "1", 2, 5 ) );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[1].minTotalCoverage, 7 );
    BOOST_REQUIRE_EQUAL( chunkedReferenceCall[1].minCoveragePerSample.size(), 1 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[1].minCoveragePerSample[0], 7 );
    BOOST_REQUIRE_EQUAL( chunkedReferenceCall[1].averageCoveragePerSample.size(), 1 );
    BOOST_CHECK_CLOSE( chunkedReferenceCall[1].averageCoveragePerSample[0], 21.0 / 3.0, 1 );
}

BOOST_AUTO_TEST_CASE( shouldComputeChunkedCoverageFor2SamplesWithAlmostConstantReadDepth )
{

    Region region = Region( "1", 0, 5 );

    wecall::io::readsummaries::coverageDeltas_t readDeltas( {{10, 0, 1, 0, 0}, {11, 1, 0, -1, 0}} );

    std::vector< readCoverage_t > chunkedReferenceCall =
        wecall::io::readsummaries::getChunkedReferenceCalls( region, 2, region.size(), readDeltas, 0.2 );

    BOOST_REQUIRE_EQUAL( chunkedReferenceCall.size(), 1 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[0].region, region );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[0].minTotalCoverage, 10 );
    BOOST_REQUIRE_EQUAL( chunkedReferenceCall[0].minCoveragePerSample.size(), 2 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[0].minCoveragePerSample[0], 10 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[0].minCoveragePerSample[1], 11 );
    BOOST_REQUIRE_EQUAL( chunkedReferenceCall[0].averageCoveragePerSample.size(), 2 );
    BOOST_CHECK_CLOSE( chunkedReferenceCall[0].averageCoveragePerSample[0], 53.0 / 5.0, 1 );
    BOOST_CHECK_CLOSE( chunkedReferenceCall[0].averageCoveragePerSample[1], 57.0 / 5.0, 1 );
}

BOOST_AUTO_TEST_CASE( shouldComputeChunkedCoverageFor2SamplesWithStrongChangesInReadDepth )
{

    Region region = Region( "1", 0, 5 );

    // expected read depth: {{5, 5, 2, 2, 2}, {6, 6, 6, 6, 0}}
    wecall::io::readsummaries::coverageDeltas_t readDeltas( {{5, 0, -3, 0, 0}, {6, 0, 0, 0, -6}} );

    std::vector< readCoverage_t > chunkedReferenceCall =
        wecall::io::readsummaries::getChunkedReferenceCalls( region, 2, region.size(), readDeltas, 0.2 );

    BOOST_REQUIRE_EQUAL( chunkedReferenceCall.size(), 3 );

    BOOST_CHECK_EQUAL( chunkedReferenceCall[0].minTotalCoverage, 5 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[0].region, Region( "1", 0, 2 ) );
    BOOST_REQUIRE_EQUAL( chunkedReferenceCall[0].minCoveragePerSample.size(), 2 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[0].minCoveragePerSample[0], 5 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[0].minCoveragePerSample[1], 6 );
    BOOST_REQUIRE_EQUAL( chunkedReferenceCall[0].averageCoveragePerSample.size(), 2 );
    BOOST_CHECK_CLOSE( chunkedReferenceCall[0].averageCoveragePerSample[0], 5.0, 1 );
    BOOST_CHECK_CLOSE( chunkedReferenceCall[0].averageCoveragePerSample[1], 6.0, 1 );

    BOOST_CHECK_EQUAL( chunkedReferenceCall[1].minTotalCoverage, 2 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[1].region, Region( "1", 2, 4 ) );
    BOOST_REQUIRE_EQUAL( chunkedReferenceCall[1].minCoveragePerSample.size(), 2 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[1].minCoveragePerSample[0], 2 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[1].minCoveragePerSample[1], 6 );
    BOOST_REQUIRE_EQUAL( chunkedReferenceCall[1].averageCoveragePerSample.size(), 2 );
    BOOST_CHECK_CLOSE( chunkedReferenceCall[1].averageCoveragePerSample[0], 2.0, 1 );
    BOOST_CHECK_CLOSE( chunkedReferenceCall[1].averageCoveragePerSample[1], 6.0, 1 );

    BOOST_CHECK_EQUAL( chunkedReferenceCall[2].minTotalCoverage, 0 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[2].region, Region( "1", 4, 5 ) );
    BOOST_REQUIRE_EQUAL( chunkedReferenceCall[2].minCoveragePerSample.size(), 2 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[2].minCoveragePerSample[0], 2 );
    BOOST_CHECK_EQUAL( chunkedReferenceCall[2].minCoveragePerSample[1], 0 );
    BOOST_REQUIRE_EQUAL( chunkedReferenceCall[2].averageCoveragePerSample.size(), 2 );
    BOOST_CHECK_CLOSE( chunkedReferenceCall[2].averageCoveragePerSample[0], 2.0, 1 );
    BOOST_CHECK_CLOSE( chunkedReferenceCall[2].averageCoveragePerSample[1], 0.0, 1 );
}