// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "utils/referenceSequence.hpp"
#include "variant/variantNormalizer.hpp"
#include <boost/optional/optional.hpp>

using wecall::utils::ReferenceSequence;
using wecall::caller::Region;
using wecall::utils::BasePairSequence;
using wecall::variant::VariantNormalizer;
using wecall::variant::Variant;
using wecall::variant::variantSet_t;
using wecall::variant::varPtr_t;

BOOST_AUTO_TEST_CASE( shouldRealWorldWrapper1 )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( "19", 7999291, 7999291 + 97 ),
        "ACGGCACCCCCAGAGCCACGAGCACACCGCCACCCCCAGAGCCACGAGCACACCGGCACCCCCAGAGCCATGAGCACACCGGCACCCCCAGAGCCAT" );

    const auto alt =
        "CCGGCACCCCCAGAGCCACGAGCACACCGGCACCCCCAGAGCCATGAGCACACCGGCACCCCCAGAGCCACGAGCACACCAGCACCCCCAGAGCCAC";

    const VariantNormalizer variantNormalizer( referenceSequence );

    const auto variants = variantNormalizer.getNormalized( referenceSequence->region(), alt, boost::none );
    BOOST_REQUIRE_EQUAL( variants.size(), 6 );
    std::vector< varPtr_t > vecVariants( variants.cbegin(), variants.cend() );

    BOOST_CHECK_EQUAL( *vecVariants[0], Variant( referenceSequence, Region( "19", 7999291, 7999292 ), "C" ) );
    BOOST_CHECK_EQUAL( *vecVariants[1], Variant( referenceSequence, Region( "19", 7999320, 7999321 ), "G" ) );
    BOOST_CHECK_EQUAL( *vecVariants[2], Variant( referenceSequence, Region( "19", 7999335, 7999336 ), "T" ) );
    BOOST_CHECK_EQUAL( *vecVariants[3], Variant( referenceSequence, Region( "19", 7999361, 7999362 ), "C" ) );
    BOOST_CHECK_EQUAL( *vecVariants[4], Variant( referenceSequence, Region( "19", 7999371, 7999372 ), "A" ) );
    BOOST_CHECK_EQUAL( *vecVariants[5], Variant( referenceSequence, Region( "19", 7999387, 7999388 ), "C" ) );
}

BOOST_AUTO_TEST_CASE( shouldRealWorldWrapper2 )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "20", 33029373, 33029373 + 36 ),
                                                                          "TACCTGCCACCATGCCCAGCTAATTTTTTTTTTTTG" );

    const auto alt = "CACCTGCCACCATGCCCAGCTAATTTTTTTTTTTTTTTTTTTG";

    const VariantNormalizer variantNormalizer( referenceSequence );

    const auto variants = variantNormalizer.getNormalized( referenceSequence->region(), alt, boost::none );
    BOOST_REQUIRE_EQUAL( variants.size(), 2 );
    std::vector< varPtr_t > vecVariants( variants.cbegin(), variants.cend() );

    BOOST_CHECK_EQUAL( *vecVariants[0], Variant( referenceSequence, Region( "20", 33029373, 33029374 ), "C" ) );
    BOOST_CHECK_EQUAL( *vecVariants[1], Variant( referenceSequence, Region( "20", 33029396, 33029396 ), "TTTTTTT" ) );
}

BOOST_AUTO_TEST_CASE( ShouldRealWorldWrapper3 )
{
    const auto contig = "22";
    const auto referenceSequence =
        std::make_shared< ReferenceSequence >( Region( contig, 35411577, 35411577 + 14 ), "AAGCCAGGTGTGGT" );

    const auto alt = "GCCAGGGGTGGG";

    const VariantNormalizer variantNormalizer( referenceSequence );

    const auto variants = variantNormalizer.getNormalized( referenceSequence->region(), alt, boost::none );
    BOOST_REQUIRE_EQUAL( variants.size(), 3 );
    std::vector< varPtr_t > vecVariants( variants.cbegin(), variants.cend() );

    BOOST_CHECK_EQUAL( *vecVariants[0], Variant( referenceSequence, Region( contig, 35411577, 35411579 ), "" ) );
    BOOST_CHECK_EQUAL( *vecVariants[1], Variant( referenceSequence, Region( contig, 35411585, 35411586 ), "G" ) );
    BOOST_CHECK_EQUAL( *vecVariants[2], Variant( referenceSequence, Region( contig, 35411590, 35411591 ), "G" ) );
}

BOOST_AUTO_TEST_CASE( ShouldRealWorldWrapper4 )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( "21", 19785531, 19785531 + 173 ),
        "TCCTCTTCCTTCCTTCCTTTTCTCTCTCTTTCTTCCTTCCTTTCTTTCTTTCTCTTTCTTTTTCTTTTTCTTTCTTTCTTTTCTTTCTTTTCTTTCTTTTCCTTCCTTCC"
        "TTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTTCCTTCCTTCCTTCCTTTCTTTCTA" );

    const auto alt =
        "TCCTCTTCCTTCCTTCCTTTTCTCTCTCTCTCTTTCTTCCTTCCTTTCTTTCTTTCTCTTTCTTTTTCTTTTTCTTTCTTTCTTTTTCTTTCTTTCTTTTCTTTCTTTTC"
        "TTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTTCCTTCCTTCCTTCCTTTCTTTCTA";

    const VariantNormalizer variantNormalizer( referenceSequence );

    const auto variants = variantNormalizer.getNormalized( referenceSequence->region(), alt, boost::none );
    BOOST_REQUIRE_EQUAL( variants.size(), 3 );
    std::vector< varPtr_t > vecVariants( variants.cbegin(), variants.cend() );

    BOOST_CHECK_EQUAL( *vecVariants[0], Variant( referenceSequence, Region( "21", 19785551, 19785551 ), "TCTC" ) );
    BOOST_CHECK_EQUAL( *vecVariants[1], Variant( referenceSequence, Region( "21", 19785610, 19785610 ), "TTTTC" ) );
    BOOST_CHECK_EQUAL( *vecVariants[2],
                       Variant( referenceSequence, Region( "21", 19785629, 19785629 ), "TTCTTTCCTTCCTTCC" ) );
}

BOOST_AUTO_TEST_CASE( ShouldRealWorldWrapper5 )
{
    const auto contig = "19";
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( contig, 46545338, 46545338 + 76 ),
        "GCCCCAGCCTCCCAAAGTGCATTGATTTTGTTGTTGTTGTGCTTATTTGCACTCCAGCCTGGCCTCTCCTTTCTTG" );

    const auto alt = "GCCCCAGCCTCCCAAAGTGCTGGGATTACAAGTGTGAACCATCGTGCCTGGCCTCTCCTTTCTTG";

    const VariantNormalizer variantNormalizer( referenceSequence );

    const auto variants = variantNormalizer.getNormalized( referenceSequence->region(), alt, boost::none );

    BOOST_REQUIRE_EQUAL( variants.size(), 9 );
    std::vector< varPtr_t > vecVariants( variants.cbegin(), variants.cend() );

    BOOST_CHECK_EQUAL( *vecVariants[0], Variant( referenceSequence, Region( contig, 46545358, 46545375 ), "" ) );
    BOOST_CHECK_EQUAL( *vecVariants[1], Variant( referenceSequence, Region( contig, 46545377, 46545378 ), "G" ) );
    BOOST_CHECK_EQUAL( *vecVariants[2], Variant( referenceSequence, Region( contig, 46545379, 46545380 ), "A" ) );
    BOOST_CHECK_EQUAL( *vecVariants[3], Variant( referenceSequence, Region( contig, 46545383, 46545383 ), "CAAG" ) );
    BOOST_CHECK_EQUAL( *vecVariants[4], Variant( referenceSequence, Region( contig, 46545384, 46545385 ), "G" ) );
    BOOST_CHECK_EQUAL( *vecVariants[5], Variant( referenceSequence, Region( contig, 46545387, 46545388 ), "A" ) );
    BOOST_CHECK_EQUAL( *vecVariants[6], Variant( referenceSequence, Region( contig, 46545390, 46545390 ), "CA" ) );
    BOOST_CHECK_EQUAL( *vecVariants[7], Variant( referenceSequence, Region( contig, 46545392, 46545393 ), "G" ) );
    BOOST_CHECK_EQUAL( *vecVariants[8], Variant( referenceSequence, Region( contig, 46545393, 46545394 ), "T" ) );
}

BOOST_AUTO_TEST_CASE( ShouldRealWorldWrapper6 )
{
    const auto contig = "22";
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( contig, 1, 1 + 42 ), "GAGGGTCCTGCAAGGAACTGCGGGAAGTCTGGAGACGGCAGG" );

    const auto alt = "AGCCCCCCACCACCCCTGCAAGGAACTGCGGGAAGTCTGGAGACGGCAGA";

    const VariantNormalizer variantNormalizer( referenceSequence );

    const auto variants = variantNormalizer.getNormalized( referenceSequence->region(), alt, boost::none );
    //    BOOST_REQUIRE_EQUAL( variants.size(), 3 );
    std::vector< varPtr_t > vecVariants( variants.cbegin(), variants.cend() );

    BOOST_CHECK_EQUAL( *vecVariants[0], Variant( referenceSequence, Region( contig, 1, 1 ), "AGCCCCCCACCACC" ) );
    BOOST_CHECK_EQUAL( *vecVariants[1], Variant( referenceSequence, Region( contig, 1, 7 ), "" ) );
    BOOST_CHECK_EQUAL( *vecVariants[2], Variant( referenceSequence, Region( contig, 42, 43 ), "A" ) );
}

BOOST_AUTO_TEST_CASE( shouldWorkWithEmptyUnnormalizedVariants )
{
    const auto contig = "owen";
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( contig, 1, 4 ), "ATG" );

    const auto alt = "ATG";

    const variantSet_t unnormalizedVariants = {};

    const VariantNormalizer variantNormalizer( referenceSequence );

    const auto variants = variantNormalizer.getNormalized( referenceSequence->region(), alt,
                                                           (boost::optional< variantSet_t >)unnormalizedVariants );
    BOOST_REQUIRE_EQUAL( variants.size(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldWorkWithNonEmptyUnnormalizedVariants )
{
    const auto contig = "owan";
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( contig, 1, 4 ), "ATG" );

    const auto alt = "AGG";

    const variantSet_t unnormalizedVariants = {
        std::make_shared< Variant >( referenceSequence, Region( contig, 2, 4 ), "" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, 2, 2 ), "G" )};

    const VariantNormalizer variantNormalizer( referenceSequence );

    const auto variants = variantNormalizer.getNormalized( referenceSequence->region(), alt,
                                                           (boost::optional< variantSet_t >)unnormalizedVariants );
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    std::vector< varPtr_t > vecVariants( variants.cbegin(), variants.cend() );

    BOOST_CHECK_EQUAL( *vecVariants[0], Variant( referenceSequence, Region( contig, 2, 3 ), "G" ) );
}

BOOST_AUTO_TEST_CASE( shouldWorkWithStrangeUnnormalizedVariants )
{
    const auto contig = "owan";
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( contig, 1, 4 ), "ATG" );

    const auto alt = "AGG";

    const variantSet_t unnormalizedVariants = {
        std::make_shared< Variant >( referenceSequence, Region( contig, 2, 4 ), "" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, 2, 2 ), "G" )};

    const VariantNormalizer variantNormalizer( referenceSequence );

    const auto variants = variantNormalizer.getNormalized( referenceSequence->region(), alt,
                                                           (boost::optional< variantSet_t >)unnormalizedVariants );
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    std::vector< varPtr_t > vecVariants( variants.cbegin(), variants.cend() );

    BOOST_CHECK_EQUAL( *vecVariants[0], Variant( referenceSequence, Region( contig, 2, 3 ), "G" ) );
}

BOOST_AUTO_TEST_CASE( shouldNotLeftAlignIntoPadding1 )
{
    const auto contig = "own";
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( contig, 1, 5 ), "AAAG" );

    const auto alt = "AAG";

    const variantSet_t unnormalizedVariants = {
        std::make_shared< Variant >( referenceSequence, Region( contig, 3, 4 ), "" )};

    const VariantNormalizer variantNormalizer( referenceSequence );

    const auto variants = variantNormalizer.getNormalized( referenceSequence->region(), alt,
                                                           (boost::optional< variantSet_t >)unnormalizedVariants );
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    std::vector< varPtr_t > vecVariants( variants.cbegin(), variants.cend() );

    BOOST_CHECK_EQUAL( *vecVariants[0], Variant( referenceSequence, Region( contig, 3, 4 ), "" ) );
}

BOOST_AUTO_TEST_CASE( shouldLeftAlignIntoPadding1 )
{
    const bool skipped = true;
    if ( skipped )
    {
        return;
    }
    const auto contig = "own";
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( contig, 1, 5 ), "AAAG" );

    const auto alt = "AAG";

    const variantSet_t unnormalizedVariants = {
        std::make_shared< Variant >( referenceSequence, Region( contig, 3, 4 ), "" )};

    const VariantNormalizer variantNormalizer( referenceSequence );

    const auto variants = variantNormalizer.getNormalized( referenceSequence->region(), alt,
                                                           (boost::optional< variantSet_t >)unnormalizedVariants );
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    std::vector< varPtr_t > vecVariants( variants.cbegin(), variants.cend() );

    BOOST_CHECK_EQUAL( *vecVariants[0], Variant( referenceSequence, Region( contig, 1, 2 ), "" ) );
}

BOOST_AUTO_TEST_CASE( shouldLeftAlignIntoPadding2 )
{
    const bool skipped = true;
    if ( skipped )
    {
        return;
    }
    const auto contig = "own";
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( contig, 1, 5 ), "AAAG" );

    const auto alt = "ATG";

    const variantSet_t unnormalizedVariants = {
        std::make_shared< Variant >( referenceSequence, Region( contig, 2, 3 ), "T" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, 3, 4 ), "" )};

    const VariantNormalizer variantNormalizer( referenceSequence );

    const auto variants = variantNormalizer.getNormalized( referenceSequence->region(), alt,
                                                           (boost::optional< variantSet_t >)unnormalizedVariants );
    BOOST_REQUIRE_EQUAL( variants.size(), 2 );
    std::vector< varPtr_t > vecVariants( variants.cbegin(), variants.cend() );

    BOOST_CHECK_EQUAL( *vecVariants[0], Variant( referenceSequence, Region( contig, 1, 2 ), "" ) );
    BOOST_CHECK_EQUAL( *vecVariants[1], Variant( referenceSequence, Region( contig, 3, 4 ), "T" ) );
}