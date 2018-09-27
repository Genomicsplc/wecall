// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "utils/NeedlemanWunsch.hpp"

using echidna::utils::BasePairSequence;
using echidna::utils::NeedlemanWunsch;
using echidna::utils::NWPenalties;
using echidna::utils::NWVariant;

BOOST_AUTO_TEST_CASE( shouldRealWorld1 )
{
    const BasePairSequence ref =
        "ACGGCACCCCCAGAGCCACGAGCACACCGCCACCCCCAGAGCCACGAGCACACCGGCACCCCCAGAGCCATGAGCACACCGGCACCCCCAGAGCCAT";
    const BasePairSequence alt =
        "CCGGCACCCCCAGAGCCACGAGCACACCGGCACCCCCAGAGCCATGAGCACACCGGCACCCCCAGAGCCACGAGCACACCAGCACCCCCAGAGCCAC";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 6 );
    BOOST_CHECK_EQUAL( nWVariants[5], NWVariant( 0, 1, "C" ) );
    BOOST_CHECK_EQUAL( nWVariants[4], NWVariant( 29, 30, "G" ) );
    BOOST_CHECK_EQUAL( nWVariants[3], NWVariant( 44, 45, "T" ) );
    BOOST_CHECK_EQUAL( nWVariants[2], NWVariant( 70, 71, "C" ) );
    BOOST_CHECK_EQUAL( nWVariants[1], NWVariant( 80, 81, "A" ) );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 96, 97, "C" ) );
}

BOOST_AUTO_TEST_CASE( shouldRealWorld2 )
{
    const BasePairSequence ref = "TACCTGCCACCATGCCCAGCTAATTTTTTTTTTTTG";
    const BasePairSequence alt = "CACCTGCCACCATGCCCAGCTAATTTTTTTTTTTTTTTTTTTG";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 2 );
    BOOST_CHECK_EQUAL( nWVariants[1], NWVariant( 0, 1, "C" ) );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 23, 23, "TTTTTTT" ) );
}

BOOST_AUTO_TEST_CASE( shouldRealWorld3 )
{
    const BasePairSequence ref = "AAGCCAGGTGTGGT";
    const BasePairSequence alt = "GCCAGGGGTGGG";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 3 );
    BOOST_CHECK_EQUAL( nWVariants[2], NWVariant( 0, 2, "" ) );
    BOOST_CHECK_EQUAL( nWVariants[1], NWVariant( 8, 9, "G" ) );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 13, 14, "G" ) );
}

BOOST_AUTO_TEST_CASE( shouldRealWorld4 )
{
    const BasePairSequence ref =
        "TCCTCTTCCTTCCTTCCTTTTCTCTCTCTTTCTTCCTTCCTTTCTTTCTTTCTCTTTCTTTTTCTTTTTCTTTCTTTCTTTTCTTTCTTTTCTTTCTTTTCCTTCCTTCC"
        "TTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTTCCTTCCTTCCTTCCTTTCTTTCTA";
    const BasePairSequence alt =
        "TCCTCTTCCTTCCTTCCTTTTCTCTCTCTCTCTTTCTTCCTTCCTTTCTTTCTTTCTCTTTCTTTTTCTTTTTCTTTCTTTCTTTTTCTTTCTTTCTTTTCTTTCTTTTC"
        "TTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTTCCTTCCTTCCTTCCTTTCTTTCTA";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 3 );
    BOOST_CHECK_EQUAL( nWVariants[2], NWVariant( 20, 20, "TCTC" ) );
    BOOST_CHECK_EQUAL( nWVariants[1], NWVariant( 79, 79, "TTTTC" ) );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 98, 98, "TTCTTTCCTTCCTTCC" ) );
}

BOOST_AUTO_TEST_CASE( shouldRealWorld5 )
{
    const BasePairSequence ref = "GCCCCAGCCTCCCAAAGTGCATTGATTTTGTTGTTGTTGTGCTTATTTGCACTCCAGCCTGGCCTCTCCTTTCTTG";
    const BasePairSequence alt = "GCCCCAGCCTCCCAAAGTGCTGGGATTACAAGTGTGAACCATCGTGCCTGGCCTCTCCTTTCTTG";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 9 );
    BOOST_CHECK_EQUAL( nWVariants[8], NWVariant( 20, 37, "" ) );
    BOOST_CHECK_EQUAL( nWVariants[7], NWVariant( 39, 40, "G" ) );
    BOOST_CHECK_EQUAL( nWVariants[6], NWVariant( 41, 42, "A" ) );
    BOOST_CHECK_EQUAL( nWVariants[5], NWVariant( 45, 45, "CAAG" ) );
    BOOST_CHECK_EQUAL( nWVariants[4], NWVariant( 46, 47, "G" ) );
    BOOST_CHECK_EQUAL( nWVariants[3], NWVariant( 49, 50, "A" ) );
    BOOST_CHECK_EQUAL( nWVariants[2], NWVariant( 52, 52, "CA" ) );
    BOOST_CHECK_EQUAL( nWVariants[1], NWVariant( 54, 55, "G" ) );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 55, 56, "T" ) );
}

BOOST_AUTO_TEST_CASE( shouldRealWorld6 )
{
    const BasePairSequence ref = "GAGGGTCCTGCAAGGAACTGCGGGAAGTCTGGAGACGGCAGG";
    const BasePairSequence alt = "AGCCCCCCACCACCCCTGCAAGGAACTGCGGGAAGTCTGGAGACGGCAGA";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 3 );
    BOOST_CHECK_EQUAL( nWVariants[2], NWVariant( 0, 0, "AGCCCCCCACCACC" ) );
    BOOST_CHECK_EQUAL( nWVariants[1], NWVariant( 0, 6, "" ) );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 41, 42, "A" ) );
}

BOOST_AUTO_TEST_CASE( shouldNoVariants )
{
    const BasePairSequence ref = "ATG";
    const BasePairSequence alt = "ATG";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldPureOneBasedSNP )
{
    const BasePairSequence ref = "A";
    const BasePairSequence alt = "T";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 0, 1, "T" ) );
}

BOOST_AUTO_TEST_CASE( shouldPureOneBasedDeletion )
{
    const BasePairSequence ref = "A";
    const BasePairSequence alt = "";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 0, 1, "" ) );
}

BOOST_AUTO_TEST_CASE( shouldPureOneBasedInsertion )
{
    const BasePairSequence ref = "";
    const BasePairSequence alt = "A";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 0, 0, "A" ) );
}

BOOST_AUTO_TEST_CASE( shouldPureMultiBasedDeletion )
{
    const BasePairSequence ref = "AAAAAAAAAA";
    const BasePairSequence alt = "";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 0, 10, "" ) );
}

BOOST_AUTO_TEST_CASE( shouldPureMultiBasedInsertion )
{
    const BasePairSequence ref = "";
    const BasePairSequence alt = "AAAAAAAAAA";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 0, 0, "AAAAAAAAAA" ) );
}

BOOST_AUTO_TEST_CASE( shouldFindSNPInCentre )
{
    const BasePairSequence ref = "ATG";
    const BasePairSequence alt = "AGG";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 1, 2, "G" ) );
}

BOOST_AUTO_TEST_CASE( shouldFindDeletionInCentre )
{
    const BasePairSequence ref = "ACGG";
    const BasePairSequence alt = "AGG";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 1, 2, "" ) );
}

BOOST_AUTO_TEST_CASE( shouldFindInsertionInCentre )
{
    const BasePairSequence ref = "AGG";
    const BasePairSequence alt = "ACGG";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 1, 1, "C" ) );
}

BOOST_AUTO_TEST_CASE( shouldFindTwoBaseDeletionInCentre )
{
    const BasePairSequence ref = "AGGC";
    const BasePairSequence alt = "AC";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 1, 3, "" ) );
}

BOOST_AUTO_TEST_CASE( shouldFindTwoBaseInsertionInCentre )
{
    const BasePairSequence ref = "AC";
    const BasePairSequence alt = "AGGC";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 1, 1, "GG" ) );
}

BOOST_AUTO_TEST_CASE( shouldFindSNPAtStart )
{
    const BasePairSequence ref = "ACG";
    const BasePairSequence alt = "TCG";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 0, 1, "T" ) );
}

BOOST_AUTO_TEST_CASE( shouldFindDeletionAtStart )
{
    const BasePairSequence ref = "ACGG";
    const BasePairSequence alt = "CGG";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 0, 1, "" ) );
}

BOOST_AUTO_TEST_CASE( shouldFindInsertionAtStart )
{
    const BasePairSequence ref = "CGG";
    const BasePairSequence alt = "ACGG";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 0, 0, "A" ) );
}

BOOST_AUTO_TEST_CASE( shouldFindMultibaseDeletionAtStart )
{
    const BasePairSequence ref = "AACGG";
    const BasePairSequence alt = "CGG";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 0, 2, "" ) );
}

BOOST_AUTO_TEST_CASE( shouldFindMultibaseInsertionAtStart )
{
    const BasePairSequence ref = "CGG";
    const BasePairSequence alt = "AACGG";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 0, 0, "AA" ) );
}

BOOST_AUTO_TEST_CASE( shouldFindSNPAtEnd )
{
    const BasePairSequence ref = "ACG";
    const BasePairSequence alt = "ACT";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 2, 3, "T" ) );
}

BOOST_AUTO_TEST_CASE( shouldFindDeletionAtEnd )
{
    const BasePairSequence ref = "ACGT";
    const BasePairSequence alt = "ACG";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 3, 4, "" ) );
}

BOOST_AUTO_TEST_CASE( shouldFindInsertionAtEnd )
{
    const BasePairSequence ref = "ACG";
    const BasePairSequence alt = "ACGT";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 3, 3, "T" ) );
}

BOOST_AUTO_TEST_CASE( shouldPreferTwoIndelsToTenSNPs )
{
    const BasePairSequence ref = "CACCATGCCCAGCTAAT";
    const BasePairSequence alt = "ACCATGCCCAGCTTAAT";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 2 );
    BOOST_CHECK_EQUAL( nWVariants[1], NWVariant( 0, 1, "" ) );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 13, 13, "T" ) );
}

BOOST_AUTO_TEST_CASE( shouldLeftAlignInsertion )
{
    const BasePairSequence ref = "ACG";
    const BasePairSequence alt = "ACCG";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 1, 1, "C" ) );
}

BOOST_AUTO_TEST_CASE( shouldLeftAlignDeletion )
{
    const BasePairSequence ref = "ACCG";
    const BasePairSequence alt = "ACG";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 1, 2, "" ) );
}

BOOST_AUTO_TEST_CASE( shouldCreateVariantsWithGapChar )
{
    const BasePairSequence ref = "ANA";
    const BasePairSequence alt = "AAA";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );
    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
    BOOST_CHECK_EQUAL( nWVariants[0], NWVariant( 1, 2, "A" ) );
}

BOOST_AUTO_TEST_CASE( shouldTimeAlgorithm100 )
{
    const bool skipped = true;
    if ( skipped )
    {
        return;
    }
    const BasePairSequence ref =
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    const BasePairSequence alt =
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );

    struct timespec start, finish;
    double elapsed;

    clock_gettime( CLOCK_MONOTONIC, &start );

    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    clock_gettime( CLOCK_MONOTONIC, &finish );

    elapsed = ( finish.tv_sec - start.tv_sec );
    elapsed += ( finish.tv_nsec - start.tv_nsec ) / 1000000000.0;
    std::cout << "Time taken for n = 10^4: " << elapsed << std::endl;

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldTimeAlgorithmMixedInsertion )
{
    const bool skipped = true;
    if ( skipped )
    {
        return;
    }
    const BasePairSequence ref =
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    const BasePairSequence alt =
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAA";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );

    struct timespec start, finish;
    double elapsed;

    clock_gettime( CLOCK_MONOTONIC, &start );

    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    clock_gettime( CLOCK_MONOTONIC, &finish );

    elapsed = ( finish.tv_sec - start.tv_sec );
    elapsed += ( finish.tv_nsec - start.tv_nsec ) / 1000000000.0;
    std::cout << "Time taken for n = 10^5: " << elapsed << std::endl;

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
}

BOOST_AUTO_TEST_CASE( shouldTimeAlgorithmMixedInsertionWithUnnormalizedVariants )
{
    const bool skipped = true;
    if ( skipped )
    {
        return;
    }
    const BasePairSequence ref =
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    const BasePairSequence alt =
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAA";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );

    struct timespec start, finish;
    double elapsed;

    clock_gettime( CLOCK_MONOTONIC, &start );

    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    clock_gettime( CLOCK_MONOTONIC, &finish );

    elapsed = ( finish.tv_sec - start.tv_sec );
    elapsed += ( finish.tv_nsec - start.tv_nsec ) / 1000000000.0;
    std::cout << "Time taken for n = 10^5: " << elapsed << std::endl;

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
}

BOOST_AUTO_TEST_CASE( shouldTimeAlgorithmMixedDeletion )
{
    const bool skipped = true;
    if ( skipped )
    {
        return;
    }
    const BasePairSequence ref =
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAA";
    const BasePairSequence alt =
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );

    struct timespec start, finish;
    double elapsed;

    clock_gettime( CLOCK_MONOTONIC, &start );

    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    clock_gettime( CLOCK_MONOTONIC, &finish );

    elapsed = ( finish.tv_sec - start.tv_sec );
    elapsed += ( finish.tv_nsec - start.tv_nsec ) / 1000000000.0;
    std::cout << "Time taken for n = 10^5: " << elapsed << std::endl;

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 1 );
}

BOOST_AUTO_TEST_CASE( shouldTimeAlgorithm1000 )
{
    const bool skipped = true;
    if ( skipped )
    {
        return;
    }
    const BasePairSequence ref =
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAA";
    const BasePairSequence alt =
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "AAAAAAAAAA";

    const auto nWPenalties = NWPenalties();
    auto needlemanWunsch = NeedlemanWunsch( ref, alt, nWPenalties );

    struct timespec start, finish;
    double elapsed;

    clock_gettime( CLOCK_MONOTONIC, &start );

    needlemanWunsch.getScoreMatrix();
    auto nWVariants = needlemanWunsch.traceBack();

    clock_gettime( CLOCK_MONOTONIC, &finish );

    elapsed = ( finish.tv_sec - start.tv_sec );
    elapsed += ( finish.tv_nsec - start.tv_nsec ) / 1000000000.0;
    std::cout << "Time taken for n = 10^6: " << elapsed << std::endl;

    BOOST_REQUIRE_EQUAL( nWVariants.size(), 0 );
}
