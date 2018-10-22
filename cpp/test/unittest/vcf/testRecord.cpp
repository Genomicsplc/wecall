// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <cstdlib>
#include <tuple>
#include <boost/test/unit_test.hpp>
#include <boost/algorithm/string.hpp>
#include "variant/type/variant.hpp"

#include "vcf/record.hpp"
#include "vcf/filterDescription.hpp"
#include "VCFTestUtils.hpp"

using wecall::variant::varPtr_t;
using wecall::variant::Variant;
using wecall::caller::Region;
using wecall::utils::ReferenceSequence;

BOOST_AUTO_TEST_CASE( should_generate_correct_SNP )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 1, 6 ), "ATCAT" );
    wecall::vcf::Record rec( "1", 2, {}, "A", {"T"}, 20.0, {}, {}, {} );
    auto variants = rec.getVariants( referenceSequence );
    BOOST_CHECK_EQUAL( variants.size(), 1 );

    auto expected_variant = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "T" );
    BOOST_CHECK( checkVariantInVector( variants, expected_variant ) );
}

BOOST_AUTO_TEST_CASE( should_generate_correct_MNP )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 1, 6 ), "ATCAT" );
    wecall::vcf::Record rec( "1", 2, {}, "ATC", {"TCG"}, 20.0, {}, {}, {} );
    auto variants = rec.getVariants( referenceSequence );
    BOOST_CHECK_EQUAL( variants.size(), 1 );

    auto expected_variant = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 4 ), "TCG" );
    BOOST_CHECK( checkVariantInVector( variants, expected_variant ) );
}

BOOST_AUTO_TEST_CASE( should_generate_correct_deletion )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 1, 6 ), "ATCAT" );
    wecall::vcf::Record rec( "1", 5, {}, "AT", {"A"}, 20.0, {}, {}, {} );
    auto variants = rec.getVariants( referenceSequence );
    BOOST_CHECK_EQUAL( variants.size(), 1 );

    auto expected_variant = std::make_shared< Variant >( referenceSequence, Region( "1", 5, 6 ), "" );
    BOOST_CHECK( checkVariantInVector( variants, expected_variant ) );
}

BOOST_AUTO_TEST_CASE( should_ignore_mal_represented_deletion )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 1, 6 ), "ATCAT" );
    wecall::vcf::Record rec( "1", 5, {}, "AT", {"T"}, 20.0, {}, {}, {} );
    auto variants = rec.getVariants( referenceSequence );
    BOOST_CHECK_EQUAL( variants.size(), 1 );

    auto expected_variant = std::make_shared< Variant >( referenceSequence, Region( "1", 4, 5 ), "" );
    BOOST_CHECK( checkVariantInVector( variants, expected_variant ) );
}

BOOST_AUTO_TEST_CASE( should_generate_correct_insertion )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 1, 6 ), "ATCAT" );
    wecall::vcf::Record rec( "1", 2, {}, "A", {"ATT"}, 20.0, {}, {}, {} );
    auto variants = rec.getVariants( referenceSequence );
    BOOST_CHECK_EQUAL( variants.size(), 1 );

    auto expected_variant = std::make_shared< Variant >( referenceSequence, Region( "1", 2, 2 ), "TT" );
    BOOST_CHECK( checkVariantInVector( variants, expected_variant ) );
}

BOOST_AUTO_TEST_CASE( test_should_read_mal_formed_insertion )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 1, 6 ), "ATCAT" );
    wecall::vcf::Record rec( "1", 5, {}, "A", {"TA"}, 20.0, {}, {}, {} );
    auto variants = rec.getVariants( referenceSequence );
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );

    auto expected_variant = std::make_shared< Variant >( referenceSequence, Region( "1", 4, 4 ), "T" );
    BOOST_CHECK( checkVariantInVector( variants, expected_variant ) );
}

BOOST_AUTO_TEST_CASE( should_generate_correct_pair_of_SNPs )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 1, 6 ), "ATCAT" );
    wecall::vcf::Record rec( "1", 2, {}, "A", {"T", "G"}, 20.0, {}, {}, {} );
    auto variants = rec.getVariants( referenceSequence );
    BOOST_CHECK_EQUAL( variants.size(), 2 );

    auto expected_variant1 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "T" );
    BOOST_CHECK( checkVariantInVector( variants, expected_variant1 ) );
    auto expected_variant2 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "G" );
    BOOST_CHECK( checkVariantInVector( variants, expected_variant2 ) );
}

BOOST_AUTO_TEST_CASE( should_ignore_complex_indel )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 1, 6 ), "ATCAT" );
    wecall::vcf::Record rec( "1", 5, {}, "AT", {"ACA"}, 20.0, {}, {}, {} );
    auto variants = rec.getVariants( referenceSequence );
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
}
