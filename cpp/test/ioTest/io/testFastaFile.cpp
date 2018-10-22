// All content Copyright (C) 2018 Genomics plc
#include "io/fastaFile.hpp"
#include "caller/region.hpp"
#include "ioFixture.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>

using namespace wecall::caller;
using namespace wecall::test;

using wecall::utils::Interval;

BOOST_FIXTURE_TEST_CASE( testFastafileContigOrder, FastaIndexFileFixture )
{
    std::vector< std::string > expectedChroms = {"1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",
                                                 "10", "11", "12", "13", "14", "15", "16", "17", "18",
                                                 "19", "20", "21", "22", "X",  "Y",  "MT"};

    auto actualChroms = fastaIndices[0]->standardContigs();
    BOOST_CHECK_EQUAL_COLLECTIONS( actualChroms.begin(), actualChroms.begin() + expectedChroms.size(),
                                   expectedChroms.begin(), expectedChroms.end() );
}

BOOST_FIXTURE_TEST_CASE( testFastafileContigLength, FastaIndexFileFixture )
{
    // Here we check that we can recover the correct chromosome lengths
    // from the FASTA index file.

    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "1" ), 249250621 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "2" ), 243199373 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "3" ), 198022430 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "4" ), 191154276 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "5" ), 180915260 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "6" ), 171115067 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "7" ), 159138663 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "8" ), 146364022 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "9" ), 141213431 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "10" ), 135534747 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "11" ), 135006516 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "12" ), 133851895 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "13" ), 115169878 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "14" ), 107349540 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "15" ), 102531392 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "16" ), 90354753 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "17" ), 81195210 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "18" ), 78077248 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "19" ), 59128983 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "20" ), 63025520 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "21" ), 48129895 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "22" ), 51304566 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "X" ), 155270560 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "Y" ), 59373566 );
    BOOST_CHECK_EQUAL( fastaIndices[0]->getContigLength( "MT" ), 16569 );
}

BOOST_FIXTURE_TEST_CASE( testFastafileContigIntervals, FastaIndexFileFixture )
{
    // Here we check that we can recover the correct chromosome lengths
    // from the FASTA index file.

    const auto contigs = fastaIndices[0]->contigs();

    BOOST_CHECK_EQUAL( contigs.at( "1" ), Interval( 0, 249250621 ) );
    BOOST_CHECK_EQUAL( contigs.at( "2" ), Interval( 0, 243199373 ) );
    BOOST_CHECK_EQUAL( contigs.at( "3" ), Interval( 0, 198022430 ) );
    BOOST_CHECK_EQUAL( contigs.at( "4" ), Interval( 0, 191154276 ) );
    BOOST_CHECK_EQUAL( contigs.at( "5" ), Interval( 0, 180915260 ) );
    BOOST_CHECK_EQUAL( contigs.at( "6" ), Interval( 0, 171115067 ) );
    BOOST_CHECK_EQUAL( contigs.at( "7" ), Interval( 0, 159138663 ) );
    BOOST_CHECK_EQUAL( contigs.at( "8" ), Interval( 0, 146364022 ) );
    BOOST_CHECK_EQUAL( contigs.at( "9" ), Interval( 0, 141213431 ) );
    BOOST_CHECK_EQUAL( contigs.at( "10" ), Interval( 0, 135534747 ) );
    BOOST_CHECK_EQUAL( contigs.at( "11" ), Interval( 0, 135006516 ) );
    BOOST_CHECK_EQUAL( contigs.at( "12" ), Interval( 0, 133851895 ) );
    BOOST_CHECK_EQUAL( contigs.at( "13" ), Interval( 0, 115169878 ) );
    BOOST_CHECK_EQUAL( contigs.at( "14" ), Interval( 0, 107349540 ) );
    BOOST_CHECK_EQUAL( contigs.at( "15" ), Interval( 0, 102531392 ) );
    BOOST_CHECK_EQUAL( contigs.at( "16" ), Interval( 0, 90354753 ) );
    BOOST_CHECK_EQUAL( contigs.at( "17" ), Interval( 0, 81195210 ) );
    BOOST_CHECK_EQUAL( contigs.at( "18" ), Interval( 0, 78077248 ) );
    BOOST_CHECK_EQUAL( contigs.at( "19" ), Interval( 0, 59128983 ) );
    BOOST_CHECK_EQUAL( contigs.at( "20" ), Interval( 0, 63025520 ) );
    BOOST_CHECK_EQUAL( contigs.at( "21" ), Interval( 0, 48129895 ) );
    BOOST_CHECK_EQUAL( contigs.at( "22" ), Interval( 0, 51304566 ) );
    BOOST_CHECK_EQUAL( contigs.at( "X" ), Interval( 0, 155270560 ) );
    BOOST_CHECK_EQUAL( contigs.at( "Y" ), Interval( 0, 59373566 ) );
    BOOST_CHECK_EQUAL( contigs.at( "MT" ), Interval( 0, 16569 ) );
}

BOOST_FIXTURE_TEST_CASE( testFastafileGetSequence, FastaFileFixture )
{
    // Here we check that we can recover some sequences from a FASTA file.

    BOOST_CHECK_EQUAL( refFiles[0]->getSequence( Region( "1", 55, 65 ) ).sequence(),
                       wecall::utils::BasePairSequence( "NNNNNCGCAG" ) );

    BOOST_CHECK_EQUAL( refFiles[0]->getSequence( Region( "2", 55, 65 ) ).sequence(),
                       wecall::utils::BasePairSequence( "NNNNNCCACA" ) );

    BOOST_CHECK_EQUAL( refFiles[0]->getSequence( Region( "1", 235, 245 ) ).sequence(),
                       wecall::utils::BasePairSequence( "CAGGANNNNN" ) );

    BOOST_CHECK_EQUAL( refFiles[0]->getSequence( Region( "2", 205, 215 ) ).sequence(),
                       wecall::utils::BasePairSequence( "ATGGANNNNN" ) );

    BOOST_CHECK_THROW( refFiles[0]->getSequence( Region( "3", 205, 215 ) ).sequence(),
                       wecall::utils::wecall_exception );
}

BOOST_FIXTURE_TEST_CASE( testFastafileGetForClusterWithPaddingWiderThanMaximumRegion, FastaFileFixture )
{
    const Region clusterRegion( "1", 110, 120 );
    const Region maximalRefFileRegion( "1", 103, 127 );
    const int64_t desiredPadding = 10;

    const auto refSequenceForCluster = refFiles[0]->getSequence( maximalRefFileRegion );
    const auto paddedSequence =
        refSequenceForCluster.getPaddedSequence( clusterRegion, maximalRefFileRegion, desiredPadding );

    BOOST_CHECK_EQUAL( paddedSequence.region(), Region( "1", 100, 130 ) );
    BOOST_CHECK_EQUAL( paddedSequence.sequence(), wecall::utils::BasePairSequence( "NNN" ) +
                                                      wecall::utils::BasePairSequence( "GTGGCGCAGGCGCAGAGAGGCGCA" ) +
                                                      wecall::utils::BasePairSequence( "NNN" ) );
}

BOOST_FIXTURE_TEST_CASE( testFastafileGetForClusterThrowsIfMaximalRegionNotContainingClusterRegion, FastaFileFixture )
{
    const auto maximumRefFileRegion = Region( "1", 111, 127 );
    const auto referenceSequence = refFiles[0]->getSequence( maximumRefFileRegion );

    BOOST_CHECK_THROW( referenceSequence.getPaddedSequence( Region( "1", 110, 120 ), maximumRefFileRegion, 10 ),
                       wecall::utils::wecall_exception );
    BOOST_CHECK_THROW( referenceSequence.getPaddedSequence( Region( "1", 111, 128 ), maximumRefFileRegion, 10 ),
                       wecall::utils::wecall_exception );
}

//-------------------------------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( testFastafileGetPaddedSequence, FastaFileFixture )
{
    // Here we check that we can recover some known sequences from
    // a FASTA file. Regions which are beyond the chromosome arAGGCGCACCG

    BOOST_CHECK_EQUAL( refFiles[0]->getSequence( Region( "1", -10, 0 ) ).sequence(),
                       wecall::utils::BasePairSequence( "NNNNNNNNNN" ) );
    BOOST_CHECK_EQUAL( refFiles[0]->getSequence( Region( "1", 240, 250 ) ).sequence(),
                       wecall::utils::BasePairSequence( "NNNNNNNNNN" ) );
}

//-------------------------------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( testFastafileCaching, FastaFileFixture )
{
    // Here we check that we can recover some known sequences from
    // a FASTA file, with caching.

    // Cache a large chunk of chr1
    refFiles[0]->cacheSequence( wecall::caller::Region( "1", 0, 240 ) );

    BOOST_CHECK_EQUAL( refFiles[0]->getSequence( Region( "1", 100, 110 ) ).sequence(),
                       wecall::utils::BasePairSequence( "GGCGTGGCGC" ) );
    BOOST_CHECK_EQUAL(
        refFiles[0]->getSequence( Region( "1", 120, 180 ) ).sequence(),
        wecall::utils::BasePairSequence( "AGGCGCACCGCGCCGGCGCAGGCGCAGAGACACATGCTAGCGCGTCCAGGGGTGGAGGCG" ) );

    // doesn't throw at the attempt to cache beyond the end of the sequence.
    refFiles[0]->cacheSequence( wecall::caller::Region( "1", 0, 20000 ) );

    // Check boundaries of cache. Can get same region as specified by cache.
    Region smallRegion( "1", 100, 110 );
    refFiles[0]->cacheSequence( smallRegion );
    BOOST_CHECK_EQUAL( refFiles[0]->getSequence( smallRegion ).sequence(),
                       wecall::utils::BasePairSequence( "GGCGTGGCGC" ) );
}
