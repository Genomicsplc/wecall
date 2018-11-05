// All content Copyright (C) 2018 Genomics plc
#include "alignment/galign.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/algorithm/string/join.hpp>
#include <string>

BOOST_AUTO_TEST_CASE( testGAlignWithHighGapOpeningPenalties3BaseReference )
{
    std::string padding = "NNNNNNNN";
    std::string haplotype = "ATGN";
    const wecall::utils::BasePairSequence haplotypeSequence = padding + haplotype + padding;
    const wecall::alignment::localGapOpenPenalties_t localGapOpen( haplotypeSequence.size(),
                                                                    std::numeric_limits< int8_t >::max() );

    wecall::utils::QualitySequence qual = {10, 10, 10};

    const short gapExtend = 1;
    const short nucleotidePrior = 4;

    const wecall::alignment::GAlign align( haplotypeSequence, gapExtend, nucleotidePrior, localGapOpen );

    BOOST_CHECK_THROW( align.computeAlignmentPhredScore( "ATG", qual, 7 ), wecall::utils::wecall_exception );

    char * aln1 = new char[30];
    char * aln2 = new char[30];

    BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( "ATG", {10, 10, 10}, 8, aln1, aln2 ), 0 );

    BOOST_CHECK_EQUAL( "ATG", aln1 );
    BOOST_CHECK_EQUAL( "ATG", aln2 );

    BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( "ATG", {10, 10, 10}, 9, aln1, aln2 ), 0 );

    BOOST_CHECK_EQUAL( "ATG", aln1 );
    BOOST_CHECK_EQUAL( "ATG", aln2 );

    BOOST_CHECK_THROW( align.computeAlignmentPhredScore( "ATG", {10, 10, 10}, 10 ), wecall::utils::wecall_exception );

    delete[] aln1;
    delete[] aln2;
}

BOOST_AUTO_TEST_CASE( testGAlignScoreIsNotAffectedByExtraMatchingBases )
{
    std::string padding = "NNNNNNNN";
    std::string haplotype = "ATAT";
    const wecall::utils::BasePairSequence haplotypeSequence = padding + haplotype + padding;
    const wecall::alignment::localGapOpenPenalties_t localGapOpen( haplotypeSequence.size(),
                                                                    std::numeric_limits< int8_t >::max() );

    const short gapExtend = 1;
    const short nucleotidePrior = 4;

    const wecall::alignment::GAlign align( haplotypeSequence, gapExtend, nucleotidePrior, localGapOpen );

    char * aln1 = new char[30];
    char * aln2 = new char[30];

    BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( "A", {10}, 8, aln1, aln2 ), 0 );

    BOOST_CHECK_EQUAL( "A", aln1 );
    BOOST_CHECK_EQUAL( "A", aln2 );

    BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( "AT", {10, 10}, 8, aln1, aln2 ), 0 );

    BOOST_CHECK_EQUAL( "AT", aln1 );
    BOOST_CHECK_EQUAL( "AT", aln2 );

    BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( "ATA", {10, 10, 10}, 8, aln1, aln2 ), 0 );

    BOOST_CHECK_EQUAL( "ATA", aln1 );
    BOOST_CHECK_EQUAL( "ATA", aln2 );

    delete[] aln1;
    delete[] aln2;
}

BOOST_AUTO_TEST_CASE( testGAlignWithHighGapOpeningPenaltiesSNPAtoT )
{
    std::string padding = "NNNNNNNN";
    std::string haplotype = "ATG";
    std::string read = "TTG";
    const wecall::utils::BasePairSequence haplotypeSequence = padding + haplotype + padding;
    const wecall::alignment::localGapOpenPenalties_t localGapOpen( haplotypeSequence.size(),
                                                                    std::numeric_limits< int8_t >::max() );

    const short gapExtend = 1;
    const short nucleotidePrior = 4;
    const short mismatchBaseThreshold = 3 * nucleotidePrior;

    const wecall::alignment::GAlign align( haplotypeSequence, gapExtend, nucleotidePrior, localGapOpen );

    char * aln1 = new char[30];
    char * aln2 = new char[30];

    for ( char mismatchBaseQuality = 0; mismatchBaseQuality < mismatchBaseThreshold; ++mismatchBaseQuality )
    {
        BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( read, {mismatchBaseQuality, 10, 10}, 8, aln1, aln2 ),
                           mismatchBaseQuality );
        BOOST_CHECK_EQUAL( haplotype, aln1 );
        BOOST_CHECK_EQUAL( read, aln2 );
    }

    BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( read, {mismatchBaseThreshold, 10, 10}, 8, aln1, aln2 ),
                       mismatchBaseThreshold );
    BOOST_CHECK_EQUAL( "NNN", aln1 );
    BOOST_CHECK_EQUAL( read, aln2 );

    BOOST_CHECK_EQUAL(
        align.computeAlignmentPhredScore( read, {std::numeric_limits< char >::max(), 10, 10}, 8, aln1, aln2 ),
        mismatchBaseThreshold );

    BOOST_CHECK_EQUAL( "NNN", aln1 );
    BOOST_CHECK_EQUAL( read, aln2 );

    delete[] aln1;
    delete[] aln2;
}

BOOST_AUTO_TEST_CASE( testGAlignWithHighGapOpeningPenaltiesSNPAtoC )
{
    std::string padding = "NNNNNNNN";
    std::string haplotype = "ATG";
    std::string read = "CTG";
    const wecall::utils::BasePairSequence haplotypeSequence = padding + haplotype + padding;
    const wecall::alignment::localGapOpenPenalties_t localGapOpen( haplotypeSequence.size(),
                                                                    std::numeric_limits< int8_t >::max() );

    const short gapExtend = 1;
    const short nucleotidePrior = 4;
    const short mismatchBaseThreshold = 3 * nucleotidePrior;

    const wecall::alignment::GAlign align( haplotypeSequence, gapExtend, nucleotidePrior, localGapOpen );

    char * aln1 = new char[30];
    char * aln2 = new char[30];

    for ( char mismatchBaseQuality = 0; mismatchBaseQuality < mismatchBaseThreshold; ++mismatchBaseQuality )
    {
        BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( read, {mismatchBaseQuality, 10, 10}, 8, aln1, aln2 ),
                           mismatchBaseQuality );
        BOOST_CHECK_EQUAL( haplotype, aln1 );
        BOOST_CHECK_EQUAL( read, aln2 );
    }

    BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( read, {mismatchBaseThreshold, 10, 10}, 8, aln1, aln2 ),
                       mismatchBaseThreshold );
    BOOST_CHECK_EQUAL( "NNN", aln1 );
    BOOST_CHECK_EQUAL( read, aln2 );

    BOOST_CHECK_EQUAL(
        align.computeAlignmentPhredScore( read, {std::numeric_limits< char >::max(), 10, 10}, 8, aln1, aln2 ),
        mismatchBaseThreshold );

    BOOST_CHECK_EQUAL( "NNN", aln1 );
    BOOST_CHECK_EQUAL( read, aln2 );

    delete[] aln1;
    delete[] aln2;
}

BOOST_AUTO_TEST_CASE( testGAlignWithHighGapOpeningPenaltiesSNPAtoG )
{
    std::string padding = "NNNNNNNN";
    std::string haplotype = "ATC";
    std::string read = "GTC";
    const wecall::utils::BasePairSequence haplotypeSequence = padding + haplotype + padding;
    const wecall::alignment::localGapOpenPenalties_t localGapOpen( haplotypeSequence.size(),
                                                                    std::numeric_limits< int8_t >::max() );

    const short gapExtend = 1;
    const short nucleotidePrior = 4;
    const short mismatchBaseThreshold = 3 * nucleotidePrior;

    const std::vector< std::size_t > goodPositions = {padding.size()};

    const wecall::alignment::GAlign align( haplotypeSequence, gapExtend, nucleotidePrior, localGapOpen );

    char * aln1 = new char[30];
    char * aln2 = new char[30];

    for ( const auto goodPosition : goodPositions )
    {
        for ( char mismatchBaseQuality = 0; mismatchBaseQuality < mismatchBaseThreshold; ++mismatchBaseQuality )
        {

            BOOST_CHECK_EQUAL(
                align.computeAlignmentPhredScore( read, {mismatchBaseQuality, 10, 10}, goodPosition, aln1, aln2 ),
                mismatchBaseQuality );
            BOOST_CHECK_EQUAL( haplotype, aln1 );
            BOOST_CHECK_EQUAL( read, aln2 );
        }

        BOOST_CHECK_EQUAL(
            align.computeAlignmentPhredScore( read, {mismatchBaseThreshold, 10, 10}, goodPosition, aln1, aln2 ),
            mismatchBaseThreshold );
        BOOST_CHECK_EQUAL( "NNN", aln1 );
        BOOST_CHECK_EQUAL( read, aln2 );

        BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( read, {std::numeric_limits< char >::max(), 10, 10},
                                                             goodPosition, aln1, aln2 ),
                           mismatchBaseThreshold );

        BOOST_CHECK_EQUAL( "NNN", aln1 );
        BOOST_CHECK_EQUAL( read, aln2 );
    }

    delete[] aln1;
    delete[] aln2;
}

BOOST_AUTO_TEST_CASE( testGAlignWithNormalGapOpeningPenaltiesSingleBaseDeletion )
{
    std::string padding = "NNNNNNNN";
    std::string haplotype = "ACTCGAGATCGTAACG";
    std::string read = "ACTCGAGTCGTAACG";
    const wecall::utils::BasePairSequence haplotypeSequence = padding + haplotype + padding;

    const wecall::alignment::errorModel_t errorModel = {
        45, 42, 41, 39, 37, 32, 28, 23, 20, 19, 17, 16, 15, 14, 13, 12, 11, 11, 10, 9, 9, 8, 8, 7, 7,
        7,  6,  6,  6,  5,  5,  5,  4,  4,  4,  3,  3,  3,  3,  2,  2,  2,  2,  2,  1, 1, 1, 1, 1};
    const auto localGapOpen = wecall::alignment::computeGapOpen( haplotypeSequence, errorModel );

    wecall::utils::QualitySequence qual = std::string( read.size(), 10 );

    const short gapExtend = 1;
    const short nucleotidePrior = 4;

    const wecall::alignment::GAlign align( haplotypeSequence, gapExtend, nucleotidePrior, localGapOpen );

    char * aln1 = new char[30];
    char * aln2 = new char[30];

    BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( read, qual, 8, aln1, aln2 ), errorModel[0] );
    BOOST_CHECK_EQUAL( haplotype, aln1 );
    BOOST_CHECK_EQUAL( "ACTCGAG-TCGTAACG", aln2 );

    delete[] aln1;
    delete[] aln2;
}

BOOST_AUTO_TEST_CASE( testGAlignWithNormalGapOpeningPenaltiesSingleBaseInsertion )
{
    std::string padding = "NNNNNNNN";
    std::string haplotype = "NACTCGAGTCGTAACG";
    std::string read = "ACTCGAGATCGTAACG";
    const wecall::utils::BasePairSequence haplotypeSequence = padding + haplotype + padding;

    const wecall::alignment::errorModel_t errorModel = {
        45, 42, 41, 39, 37, 32, 28, 23, 20, 19, 17, 16, 15, 14, 13, 12, 11, 11, 10, 9, 9, 8, 8, 7, 7,
        7,  6,  6,  6,  5,  5,  5,  4,  4,  4,  3,  3,  3,  3,  2,  2,  2,  2,  2,  1, 1, 1, 1, 1};
    const auto localGapOpen = wecall::alignment::computeGapOpen( haplotypeSequence, errorModel );

    wecall::utils::QualitySequence qual = std::string( read.size(), 10 );

    const short gapExtend = 1;
    const short nucleotidePrior = 4;

    const wecall::alignment::GAlign align( haplotypeSequence, gapExtend, nucleotidePrior, localGapOpen );

    char * aln1 = new char[30];
    char * aln2 = new char[30];

    BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( read, qual, 8, aln1, aln2 ), errorModel[0] + nucleotidePrior );
    BOOST_CHECK_EQUAL( "ACTCGAG-TCGTAACG", aln1 );
    BOOST_CHECK_EQUAL( read, aln2 );

    delete[] aln1;
    delete[] aln2;
}

BOOST_AUTO_TEST_CASE( testGAlignWithNormalGapOpeningPenaltiesMultibaseDeletion )
{
    std::string padding = "NNNNNNNN";
    std::string haplotype = "ATACTCGAGAGCGTATCGTAACG";
    std::string read = "ATACTCGAGTCGTAACG";
    const wecall::utils::BasePairSequence haplotypeSequence = padding + haplotype + padding;

    const wecall::alignment::errorModel_t errorModel = {
        45, 42, 41, 39, 37, 32, 28, 23, 20, 19, 17, 16, 15, 14, 13, 12, 11, 11, 10, 9, 9, 8, 8, 7, 7,
        7,  6,  6,  6,  5,  5,  5,  4,  4,  4,  3,  3,  3,  3,  2,  2,  2,  2,  2,  1, 1, 1, 1, 1};
    const auto localGapOpen = wecall::alignment::computeGapOpen( haplotypeSequence, errorModel );

    wecall::utils::QualitySequence qual = std::string( read.size(), 10 );

    const short gapExtend = 1;
    const short nucleotidePrior = 4;

    const wecall::alignment::GAlign align( haplotypeSequence, gapExtend, nucleotidePrior, localGapOpen );

    char * aln1 = new char[30];
    char * aln2 = new char[30];
    const short gapLength = 6;

    BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( read, qual, 8, aln1, aln2 ),
                       errorModel[0] + ( gapLength - 1 ) * gapExtend );
    BOOST_CHECK_EQUAL( haplotype, aln1 );
    BOOST_CHECK_EQUAL( "ATACTCGAG------TCGTAACG", aln2 );

    delete[] aln1;
    delete[] aln2;
}

BOOST_AUTO_TEST_CASE( testGAlignWithNormalGapOpeningPenaltiesMultibaseInsertion )
{
    std::string padding = "NNNNNNNN";
    std::string haplotype = "ATACTCGAGTCGTAACGATACTG";
    std::string read = "ATACTCGAGAGCGTATCGTAACG";
    const wecall::utils::BasePairSequence haplotypeSequence = padding + haplotype + padding;

    const wecall::alignment::errorModel_t errorModel = {
        45, 42, 41, 39, 37, 32, 28, 23, 20, 19, 17, 16, 15, 14, 13, 12, 11, 11, 10, 9, 9, 8, 8, 7, 7,
        7,  6,  6,  6,  5,  5,  5,  4,  4,  4,  3,  3,  3,  3,  2,  2,  2,  2,  2,  1, 1, 1, 1, 1};
    const auto localGapOpen = wecall::alignment::computeGapOpen( haplotypeSequence, errorModel );

    wecall::utils::QualitySequence qual = std::string( read.size(), 10 );

    const short gapExtend = 1;
    const short nucleotidePrior = 4;

    const wecall::alignment::GAlign align( haplotypeSequence, gapExtend, nucleotidePrior, localGapOpen );

    char * aln1 = new char[30];
    char * aln2 = new char[30];
    const short gapLength = 6;

    BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( read, qual, 8, aln1, aln2 ),
                       errorModel[0] + ( gapLength - 1 ) * gapExtend + gapLength * nucleotidePrior );
    BOOST_CHECK_EQUAL( "ATACTCGAG------TCGTAACG", aln1 );
    BOOST_CHECK_EQUAL( read, aln2 );

    delete[] aln1;
    delete[] aln2;
}

BOOST_AUTO_TEST_CASE( testComputeGapOpenSymmetry )
{
    wecall::utils::BasePairSequence basePairSequence( "TTAAATAAAAAAATAAAAAAATAAATT" );

    const wecall::alignment::errorModel_t errorModel = {
        45, 42, 41, 39, 37, 32, 28, 23, 20, 19, 17, 16, 15, 14, 13, 12, 11, 11, 10, 9, 9, 8, 8, 7, 7,
        7,  6,  6,  6,  5,  5,  5,  4,  4,  4,  3,  3,  3,  3,  2,  2,  2,  2,  2,  1, 1, 1, 1, 1};

    const auto result = wecall::alignment::computeGapOpen( basePairSequence, errorModel );

    const wecall::alignment::localGapOpenPenalties_t referenceResult = {
        42, 45, 41, 42, 45, 45, 28, 32, 37, 39, 41, 42, 45, 45, 28, 32, 37, 39, 41, 42, 45, 45, 41, 42, 45, 42, 45};

    BOOST_REQUIRE_EQUAL( referenceResult.size(), result.size() );
    for ( size_t i = 0; i < referenceResult.size(); ++i )
    {
        BOOST_REQUIRE_EQUAL( referenceResult[i], result[i] );
    }
}

BOOST_AUTO_TEST_CASE( testGAlignWithVaryingPos )
{
    std::string padding = "NNNNNNNN";
    std::string noG = "AAAAAAAAAA";
    std::string haplotype = noG + "G" + noG;
    const wecall::utils::BasePairSequence haplotypeSequence = padding + haplotype + padding;
    const wecall::alignment::localGapOpenPenalties_t localGapOpen( haplotypeSequence.size(),
                                                                    std::numeric_limits< int8_t >::max() );

    wecall::utils::QualitySequence qual = {10};

    const short gapExtend = 1;
    const short nucleotidePrior = 4;

    const wecall::alignment::GAlign align( haplotypeSequence, gapExtend, nucleotidePrior, localGapOpen );

    const int paddingLength = padding.size();
    const int noGLength = noG.size();

    BOOST_CHECK_THROW( align.computeAlignmentPhredScore( "G", qual, paddingLength - 1 ),
                       wecall::utils::wecall_exception );

    char * aln1 = new char[30];
    char * aln2 = new char[30];

    for ( auto startpos = paddingLength; startpos <= noGLength; ++startpos )
    {
        BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( "G", qual, startpos, aln1, aln2 ), nucleotidePrior );

        BOOST_CHECK_EQUAL( "N", aln1 );
        BOOST_CHECK_EQUAL( "G", aln2 );
    }

    for ( auto startpos = noGLength + 1; startpos <= noGLength + 2 * paddingLength; ++startpos )
    {
        BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( "G", qual, startpos, aln1, aln2 ), 0 );

        BOOST_CHECK_EQUAL( "G", aln1 );
        BOOST_CHECK_EQUAL( "G", aln2 );
    }

    for ( auto startpos = noGLength + 2 * paddingLength + 1; startpos <= 2 * noGLength + paddingLength; ++startpos )
    {
        BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( "G", qual, startpos, aln1, aln2 ), nucleotidePrior );

        BOOST_CHECK_EQUAL( "N", aln1 );
        BOOST_CHECK_EQUAL( "G", aln2 );
    }

    BOOST_CHECK_THROW( align.computeAlignmentPhredScore( "G", qual, 2 * noGLength + paddingLength + 1 ),
                       wecall::utils::wecall_exception );

    delete[] aln1;
    delete[] aln2;
}

BOOST_AUTO_TEST_CASE( testGAlignNucleotidePrior )
{
    std::string padding = "NNNNNNNN";
    std::string haplotype = "N";
    std::string read = "A";
    const wecall::utils::BasePairSequence haplotypeSequence = padding + haplotype + padding;

    const wecall::alignment::errorModel_t errorModel = {
        45, 42, 41, 39, 37, 32, 28, 23, 20, 19, 17, 16, 15, 14, 13, 12, 11, 11, 10, 9, 9, 8, 8, 7, 7,
        7,  6,  6,  6,  5,  5,  5,  4,  4,  4,  3,  3,  3,  3,  2,  2,  2,  2,  2,  1, 1, 1, 1, 1};
    const auto localGapOpen = wecall::alignment::computeGapOpen( haplotypeSequence, errorModel );

    wecall::utils::QualitySequence qual = {10};

    const short gapExtend = 1;
    const short nucleotidePrior = 4;

    const wecall::alignment::GAlign align( haplotypeSequence, gapExtend, nucleotidePrior, localGapOpen );

    char * aln1 = new char[30];
    char * aln2 = new char[30];

    BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( read, qual, 8, aln1, aln2 ), nucleotidePrior );
    BOOST_CHECK_EQUAL( haplotype, aln1 );
    BOOST_CHECK_EQUAL( read, aln2 );

    delete[] aln1;
    delete[] aln2;
}

BOOST_AUTO_TEST_CASE( testGAlignWithEmptySequence )
{
    std::string padding = "NNNNNNNN";
    const wecall::utils::BasePairSequence haplotypeSequence = padding + "ATG" + padding;
    const wecall::alignment::localGapOpenPenalties_t localGapOpen( haplotypeSequence.size(),
                                                                    std::numeric_limits< int8_t >::max() );

    wecall::utils::QualitySequence qual = {};

    const short gapExtend = 1;
    const short nucleotidePrior = 4;

    const wecall::alignment::GAlign align( haplotypeSequence, gapExtend, nucleotidePrior, localGapOpen );

    BOOST_CHECK_EQUAL( align.computeAlignmentPhredScore( "", qual, 8 ), 0 );
}

BOOST_AUTO_TEST_CASE( testGAlignWithWrongQualLength )
{
    std::string padding = "NNNNNNNN";
    std::string haplotype = "ATG";
    std::string read = "ATG";
    const wecall::utils::BasePairSequence haplotypeSequence = padding + haplotype + padding;
    const wecall::alignment::localGapOpenPenalties_t localGapOpen( haplotypeSequence.size(),
                                                                    std::numeric_limits< int8_t >::max() );

    wecall::utils::QualitySequence qual = std::string( read.size() - 1, 10 );
    ;

    const short gapExtend = 1;
    const short nucleotidePrior = 4;

    const wecall::alignment::GAlign align( haplotypeSequence, gapExtend, nucleotidePrior, localGapOpen );

    char * aln1 = new char[30];
    char * aln2 = new char[30];

    BOOST_CHECK_THROW( align.computeAlignmentPhredScore( read, qual, 8, aln1, aln2 ),
                       wecall::utils::wecall_exception );

    delete[] aln1;
    delete[] aln2;
}
