// All content Copyright (C) 2018 Genomics plc
#include "io/fastaFile.hpp"
#include "readrecalibration/kmerDistribution.hpp"
#include "readrecalibration/commonTypes.hpp"
#include "alignment/cigarItems.hpp"
#include "caller/region.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "ioTest/io/ioFixture.hpp"

using wecall::alignment::Cigar;

BOOST_FIXTURE_TEST_CASE( testConstructionOfReferenceKmers, wecall::test::FastaFileFixture )
{
    std::string chrom( "1" );
    const int readStart = 60;
    const int readEnd = 120;

    // Given the kmerDistrubution build with this sequence,
    const auto expectedReference =
        wecall::utils::BasePairSequence( "N" ) +
        wecall::utils::BasePairSequence( "CGCAGGCGCAGAGACACATGCTACCGCGTCCAGGGGTGGAGGCGTGGCGCAGGCGCAGAG" ) +
        wecall::utils::BasePairSequence( "AGGCGCAC" );
    BOOST_CHECK_EQUAL( refFiles[0]
                           ->getSequence( wecall::caller::Region(
                               chrom, readStart - wecall::corrector::padding,
                               readEnd + wecall::corrector::padding + wecall::corrector::kmerSize ) )
                           .sequence(),
                       expectedReference );

    wecall::corrector::KmerDistribution kmerDistribution( chrom, *refFiles[0], readStart, readEnd );

    for ( int pos = 0; pos < readEnd - readStart; ++pos )
    {
        auto siteDistrubution = kmerDistribution.getSiteKmerDistribution( readStart + pos );

        auto referenceKmer = siteDistrubution.getReferenceKmer();

        auto referenceKmerAsString = wecall::corrector::show_string( referenceKmer );
        auto expectedKmer =
            expectedReference.str().substr( wecall::corrector::padding + pos, wecall::corrector::kmerSize );

        BOOST_CHECK_EQUAL( referenceKmerAsString, expectedKmer );
    }
}

BOOST_FIXTURE_TEST_CASE( testConstructionOfExtendedReferenceKmers, wecall::test::FastaFileFixture )
{
    std::string chrom( "1" );
    const int readStart = 60;
    const int readEnd = 120;

    // Given the kmerDistrubution build with this sequence,
    const auto expectedReference =
        wecall::utils::BasePairSequence( "N" ) +
        wecall::utils::BasePairSequence( "CGCAGGCGCAGAGACACATGCTACCGCGTCCAGGGGTGGAGGCGTGGCGCAGGCGCAGAG" ) +
        wecall::utils::BasePairSequence( "AGGCGCAC" );

    BOOST_CHECK_EQUAL( refFiles[0]
                           ->getSequence( wecall::caller::Region(
                               chrom, readStart - wecall::corrector::padding,
                               readEnd + wecall::corrector::padding + wecall::corrector::kmerSize ) )
                           .sequence(),
                       expectedReference );

    wecall::corrector::KmerDistribution kmerDistribution( chrom, *refFiles[0], readStart, readEnd );

    for ( int pos = 0; pos < readEnd - readStart; ++pos )
    {
        auto siteDistrubution = kmerDistribution.getSiteKmerDistribution( readStart + pos );

        auto extendedReferenceKmer = siteDistrubution.getExtReferenceKmer();

        auto extendedReferenceKmerAsString = wecall::corrector::show_string( extendedReferenceKmer );
        auto expectedKmer =
            expectedReference.str().substr( pos, wecall::corrector::kmerSize + 2 * wecall::corrector::padding );

        BOOST_CHECK_EQUAL( extendedReferenceKmerAsString, expectedKmer );
    }
}

BOOST_FIXTURE_TEST_CASE( testUpdatingHistogramWithRead, wecall::test::FastaFileFixture )
{
    std::string chrom( "1" );
    const int readStart = 60;
    const int readEnd = 70;

    // Given the kmerDistrubution build with this sequence,
    const auto expectedReference = wecall::utils::BasePairSequence( "N" ) +
                                   wecall::utils::BasePairSequence( "CGCAGGCGCA" ) +
                                   wecall::utils::BasePairSequence( "GAGACACA" );

    auto referenceSequence = refFiles[0]->getSequence(
        wecall::caller::Region( chrom, readStart - wecall::corrector::padding,
                                 readEnd + wecall::corrector::padding + wecall::corrector::kmerSize ) );
    BOOST_CHECK_EQUAL( referenceSequence.sequence(), expectedReference );

    wecall::corrector::KmerDistribution kmerDistribution( chrom, *refFiles[0], readStart, readEnd );

    for ( int i = readStart; i < readEnd; ++i )
    {
        BOOST_CHECK_EQUAL( kmerDistribution.getSiteKmerDistribution( i ).size(), 0 );
    }

    wecall::io::readPtr_t readPtr = std::make_shared< wecall::io::Read >(
        "TACCACACGT", "WHYNOTTEST", "test", Cigar( std::to_string( readEnd - readStart ) + "M" ), 0, readStart,
        BAM_FPROPER_PAIR, 40, 2000500, 200, 0,
        std::make_shared< wecall::utils::ReferenceSequence >( referenceSequence ) );

    kmerDistribution.updateKmerHistogram( readPtr );

    for ( std::size_t i = readStart; i < readEnd - wecall::corrector::kmerSize + 1; ++i )
    {
        BOOST_CHECK_EQUAL( kmerDistribution.getSiteKmerDistribution( i ).size(), 1 );
        for ( auto kmerCountPair : kmerDistribution.getSiteKmerDistribution( i ).kmerCount() )
        {
            BOOST_CHECK_EQUAL( kmerCountPair.second, 1 );
            auto expectedSequence = readPtr->sequence().substr( i - readStart, wecall::corrector::kmerSize );
            auto kmerSequence = wecall::corrector::show_string( kmerCountPair.first );
            BOOST_CHECK_EQUAL( kmerSequence, expectedSequence.str() );
        }
    }
}

BOOST_FIXTURE_TEST_CASE( testFinalizeWithNoReads, wecall::test::FastaFileFixture )
{
    std::string chrom( "1" );
    const int readStart = 60;
    const int readEnd = 120;
    const double pref = 0.05;

    // Given the kmerDistrubution build with this sequence,
    const auto expectedReference =
        wecall::utils::BasePairSequence( "N" ) +
        wecall::utils::BasePairSequence( "CGCAGGCGCAGAGACACATGCTACCGCGTCCAGGGGTGGAGGCGTGGCGCAGGCGCAGAG" ) +
        wecall::utils::BasePairSequence( "AGGCGCAC" );

    BOOST_CHECK_EQUAL( refFiles[0]
                           ->getSequence( wecall::caller::Region(
                               chrom, readStart - wecall::corrector::padding,
                               readEnd + wecall::corrector::padding + wecall::corrector::kmerSize ) )
                           .sequence(),
                       expectedReference );

    wecall::corrector::KmerDistribution kmerDistribution( chrom, *refFiles[0], readStart, readEnd );
    kmerDistribution.finalise( pref );

    for ( auto pos = kmerDistribution.start(); pos != kmerDistribution.end(); ++pos )
    {
        const auto & siteDistrubution = kmerDistribution.getSiteKmerDistribution( pos );

        const auto & priors = siteDistrubution.kmerPrior();

        BOOST_CHECK_EQUAL( priors.size(), 1 );

        for ( const auto & priorPair : priors )
        {
            BOOST_CHECK_EQUAL( wecall::corrector::show_string( priorPair.first ),
                               wecall::corrector::show_string( siteDistrubution.getReferenceKmer() ) );
            BOOST_CHECK_CLOSE( priorPair.second, 1.0, 1e-8 );
        }
    }
}

BOOST_FIXTURE_TEST_CASE( testFinalizeWithThreeMatchingReads, wecall::test::FastaFileFixture )
{
    // Test should fail if computation priors is altered,
    std::string chrom( "1" );
    const int readStart = 60;
    const int readEnd = 70;
    const int refStart = 60;
    const int refEnd = 120;
    const double pref = 0.95;

    wecall::corrector::KmerDistribution kmerDistribution( chrom, *refFiles[0], refStart, refEnd );

    auto refSequence = refFiles[0]->getSequence( wecall::caller::Region( chrom, readStart, readEnd ) );
    BOOST_CHECK_EQUAL( refSequence.sequence().str(), "CGCAGGCGCA" );
    wecall::io::readPtr_t readPtr = std::make_shared< wecall::io::Read >(
        "TGCAGGCGCA", "zzzzzzzzzz", "test", Cigar( "10M" ), 0, readStart, BAM_FPROPER_PAIR, 40, 500, 200, 0,
        std::make_shared< wecall::utils::ReferenceSequence >( refSequence ) );

    // Some reason adding more reads doesn't change the prior for the kmer with the snp. GL: This is expected.
    auto readCover = 2;
    for ( auto i = 0; i < readCover; ++i )
    {
        kmerDistribution.updateKmerHistogram( readPtr );
    }

    kmerDistribution.finalise( pref );

    const auto & siteDistrubution = kmerDistribution.getSiteKmerDistribution( readStart );
    const auto & priors = siteDistrubution.kmerPrior();
    BOOST_CHECK_EQUAL( priors.size(), 2 );
    BOOST_CHECK_EQUAL( wecall::corrector::show_string( priors[0].first ),
                       wecall::corrector::show_string( siteDistrubution.getReferenceKmer() ) );
    BOOST_CHECK_CLOSE( priors[0].second, 0.95238095238095233, 1e-5 );

    BOOST_CHECK_EQUAL(
        wecall::corrector::show_string( priors[1].first ),
        wecall::corrector::show_string( readPtr->sequence().substr( 0, wecall::corrector::kmerSize ) ) );
    BOOST_CHECK_CLOSE( priors[1].second, 0.047619047619047658, 1e-5 );

    BOOST_CHECK_CLOSE( priors[0].second + priors[1].second, 1.0, 1e-8 );
}
