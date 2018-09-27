// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <string>

#include "variant/genotype.hpp"
#include "variant/haplotype.hpp"
#include "variant/type/variant.hpp"
#include "caller/diploid/diploidAnnotate.hpp"

using echidna::variant::GenotypeVector;
using echidna::variant::HaplotypeVector;
using echidna::variant::genotypePtr_t;
using echidna::variant::varPtr_t;
using echidna::variant::variantSet_t;
using echidna::variant::Haplotype;
using echidna::utils::ReferenceSequence;
using echidna::caller::Region;
using echidna::variant::Variant;

BOOST_AUTO_TEST_CASE( testConstructorThrowsIfProvidedWithUnMergedHaplotypes )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "T" );

    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );

    haplotypeVector.push_back( variantSet_t{} );
    haplotypeVector.push_back( variantSet_t{} );

    BOOST_CHECK_THROW( GenotypeVector( 2, haplotypeVector, {0, 1}, {} ), echidna::utils::echidna_exception );
}

BOOST_AUTO_TEST_CASE( testConstructsZeroGenotypesForEmptyHaplotypeVector )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "T" );

    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );

    const GenotypeVector genotypeVector( 2, haplotypeVector, {}, {} );

    BOOST_CHECK_EQUAL( genotypeVector.size(), 0 );
}

BOOST_AUTO_TEST_CASE( testConstructsOneGenotypeForHaplotypeVectorWithOneEmptyElement )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "T" );

    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );
    haplotypeVector.push_back( variantSet_t{} );

    const auto ploidy = 2;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0}, {} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 1 );

    BOOST_CHECK_EQUAL( genotypeVector.ploidy(), ploidy );

    const auto haplotypeIndicies = genotypeVector.getHaplotypeIndices( 0 );
    BOOST_REQUIRE_EQUAL( haplotypeIndicies.size(), 2 );
    BOOST_CHECK_EQUAL( haplotypeIndicies[0], 0 );
    BOOST_CHECK_EQUAL( haplotypeIndicies[1], 0 );
}

BOOST_AUTO_TEST_CASE( testPhaseQualityComputationEqualLikelihoods )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );

    varPtr_t snp = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );
    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );
    haplotypeVector.push_back( variantSet_t{} );
    haplotypeVector.push_back( variantSet_t{snp} );

    const auto ploidy = 2;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1}, {snp} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 3 );

    // Genotype quality should be 1.0 - max-likelihood / total-likelihoods.
    const std::vector< double > genotypeLikelihoods = {0.1, 0.1, 0.1};

    // Middle genotype likelihood has two possible combinations.
    const double expectedPhaseQuality = 0;  // 1.0 - ( 0.1 * 2 ) / ( 0.1 * 2 );
    const std::size_t calledGenotypeIndex = 1;

    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( echidna::caller::annotate::computePhaseQuality(
                           calledGenotypeIndex, genotypeVector, genotypeLikelihoods ) ),
                       echidna::stats::roundPhred( echidna::stats::toPhredQ( expectedPhaseQuality ) ) );
}

BOOST_AUTO_TEST_CASE( testPhaseQualityComputationUnequalLikelihoods )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );

    varPtr_t snp = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );
    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );
    haplotypeVector.push_back( variantSet_t{} );
    haplotypeVector.push_back( variantSet_t{snp} );

    const auto ploidy = 2;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1}, {snp} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 3 );

    const std::vector< double > genotypeLikelihoods = {0.9, 0.1, 0.1};

    // Middle genotype likelihood has two possible combinations.
    const double expectedPhaseQuality = 0;  // 1.0 - 0.9 / ( 0.9 );
    const std::size_t calledGenotypeIndex = 0;

    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( echidna::caller::annotate::computePhaseQuality(
                           calledGenotypeIndex, genotypeVector, genotypeLikelihoods ) ),
                       echidna::stats::roundPhred( echidna::stats::toPhredQ( expectedPhaseQuality ) ) );
}

BOOST_AUTO_TEST_CASE( testPhaseQualityComputationSmallLikelihoods )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );

    varPtr_t snp = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );
    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );
    haplotypeVector.push_back( variantSet_t{} );
    haplotypeVector.push_back( variantSet_t{snp} );

    const auto ploidy = 2;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1}, {snp} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 3 );

    const std::vector< double > genotypeLikelihoods = {3e-9, 1e-9, 1e-9};

    // Middle genotype likelihood has two possible combinations.
    const double expectedPhaseQuality = 1.0 - 3.0 / ( 3.0 );
    const std::size_t calledGenotypeIndex = 0;

    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( echidna::caller::annotate::computePhaseQuality(
                           calledGenotypeIndex, genotypeVector, genotypeLikelihoods ) ),
                       echidna::stats::roundPhred( echidna::stats::toPhredQ( expectedPhaseQuality ) ) );
}

BOOST_AUTO_TEST_CASE( testPhaseQualityComputationLargeLikelihoods )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );

    varPtr_t snp = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );
    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );
    haplotypeVector.push_back( variantSet_t{} );
    haplotypeVector.push_back( variantSet_t{snp} );

    const auto ploidy = 2;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1}, {snp} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 3 );

    const std::vector< double > genotypeLikelihoods = {3e9, 1e9, 1e9};

    // Middle genotype likelihood has two possible combinations.
    const double expectedPhaseQuality = 1.0 - 3.0 / ( 3.0 );
    const std::size_t calledGenotypeIndex = 0;

    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( echidna::caller::annotate::computePhaseQuality(
                           calledGenotypeIndex, genotypeVector, genotypeLikelihoods ) ),
                       echidna::stats::roundPhred( echidna::stats::toPhredQ( expectedPhaseQuality ) ) );
}

BOOST_AUTO_TEST_CASE( testGenotypeQualityComputationEqualLikelihoods )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );

    varPtr_t snp = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );
    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );
    haplotypeVector.push_back( variantSet_t{} );
    haplotypeVector.push_back( variantSet_t{snp} );

    const auto ploidy = 2;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1}, {} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 3 );

    // Genotype quality should be 1.0 - max-likelihood / total-likelihoods.
    const std::vector< double > genotypeLikelihoods = {0.1, 0.1, 0.1};

    // Middle genotype likelihood has two possible combinations.
    const double expectedGenotypeQuality = 1.0 - ( 0.1 * 2 ) / ( 0.1 + 0.1 * 2 + 0.1 );

    BOOST_CHECK_EQUAL( echidna::stats::roundPhred(
                           echidna::caller::annotate::computeGenotypeQuality( genotypeVector, genotypeLikelihoods ) ),
                       echidna::stats::roundPhred( echidna::stats::toPhredQ( expectedGenotypeQuality ) ) );
}

BOOST_AUTO_TEST_CASE( testGenotypeQualityComputationUnequalLikelihoods )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );

    varPtr_t snp = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );
    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );
    haplotypeVector.push_back( variantSet_t{} );
    haplotypeVector.push_back( variantSet_t{snp} );

    const auto ploidy = 2;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1}, {snp} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 3 );

    const std::vector< double > genotypeLikelihoods = {0.9, 0.1, 0.1};

    // Middle genotype likelihood has two possible combinations.
    const double expectedGenotypeQuality = 1.0 - 0.9 / ( 0.9 + 0.1 * 2 + 0.1 );

    BOOST_CHECK_EQUAL( echidna::stats::roundPhred(
                           echidna::caller::annotate::computeGenotypeQuality( genotypeVector, genotypeLikelihoods ) ),
                       echidna::stats::roundPhred( echidna::stats::toPhredQ( expectedGenotypeQuality ) ) );
}

BOOST_AUTO_TEST_CASE( testGenotypeQualityComputationSmallLikelihoods )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );

    varPtr_t snp = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );
    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );
    haplotypeVector.push_back( variantSet_t{} );
    haplotypeVector.push_back( variantSet_t{snp} );

    const auto ploidy = 2;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1}, {} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 3 );

    const std::vector< double > genotypeLikelihoods = {3e-9, 1e-9, 1e-9};

    // Middle genotype likelihood has two possible combinations.
    const double expectedGenotypeQuality = 1.0 - 3.0 / ( 3.0 + 1.0 * 2 + 1.0 );

    BOOST_CHECK_EQUAL( echidna::stats::roundPhred(
                           echidna::caller::annotate::computeGenotypeQuality( genotypeVector, genotypeLikelihoods ) ),
                       echidna::stats::roundPhred( echidna::stats::toPhredQ( expectedGenotypeQuality ) ) );
}

BOOST_AUTO_TEST_CASE( testGenotypeQualityComputationLargeLikelihoods )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );

    varPtr_t snp = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );
    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );
    haplotypeVector.push_back( variantSet_t{} );
    haplotypeVector.push_back( variantSet_t{snp} );

    const auto ploidy = 2;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1}, {} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 3 );

    const std::vector< double > genotypeLikelihoods = {3e9, 1e9, 1e9};

    // Middle genotype likelihood has two possible combinations.
    const double expectedGenotypeQuality = 1.0 - 3.0 / ( 3.0 + 1.0 * 2 + 1.0 );

    BOOST_CHECK_EQUAL( echidna::stats::roundPhred(
                           echidna::caller::annotate::computeGenotypeQuality( genotypeVector, genotypeLikelihoods ) ),
                       echidna::stats::roundPhred( echidna::stats::toPhredQ( expectedGenotypeQuality ) ) );
}

BOOST_AUTO_TEST_CASE( testGenotypeLikelihoodsComputationEqualLikelihoods )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );

    varPtr_t snp = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );
    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );
    haplotypeVector.push_back( variantSet_t{} );
    haplotypeVector.push_back( variantSet_t{snp} );

    const auto ploidy = 2;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1}, {} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 3 );

    // Genotype quality should be 1.0 - max-likelihood / total-likelihoods.
    const std::vector< double > genotypeLikelihoods = {0.1, 0.1, 0.1};
    const auto likelihoods = echidna::caller::annotate::get_RR_RA_AA_Likelihoods_as_phred_scores(
        std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" ), haplotypeVector, genotypeVector,
        genotypeLikelihoods );

    BOOST_REQUIRE_EQUAL( likelihoods.size(), 3 );  // In ploidy = 2 expect 2 + 1 likelihoods.

    // Expect 0.1, 0.2, 0.1 from raw likelihoods as het call should be scaled by 2.0
    // This is equivalent to 0.5, 1.0, 0.5 which in phred space is rounded to 3, 0, 3
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[0] ), 3 );
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[1] ), 0 );
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[2] ), 3 );
}

BOOST_AUTO_TEST_CASE( testGenotypeLikelihoodsComputationUnequalLikelihoods )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );

    varPtr_t snp = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );
    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );
    haplotypeVector.push_back( variantSet_t{} );
    haplotypeVector.push_back( variantSet_t{snp} );

    const auto ploidy = 2;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1}, {} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 3 );

    const std::vector< double > genotypeLikelihoods = {0.9, 0.1, 0.1};
    const auto likelihoods = echidna::caller::annotate::get_RR_RA_AA_Likelihoods_as_phred_scores(
        std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" ), haplotypeVector, genotypeVector,
        genotypeLikelihoods );

    BOOST_REQUIRE_EQUAL( likelihoods.size(), 3 );  // In ploidy = 2 expect 2 + 1 likelihoods.

    // Expect 0.9, 0.2, 0.1 from raw likelihoods as het call should be scaled by 2.0
    // Equivalent to 1.0, 0.22222, 0.11111 --> 0, 6.5321, 9.54245
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[0] ), 0 );
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[1] ), 7 );
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[2] ), 10 );
}

BOOST_AUTO_TEST_CASE( testGenotypeLikelihoodsComputationSmallLikelihoods )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );

    varPtr_t snp = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );
    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );
    haplotypeVector.push_back( variantSet_t{} );
    haplotypeVector.push_back( variantSet_t{snp} );

    const auto ploidy = 2;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1}, {} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 3 );
    const std::vector< double > genotypeLikelihoods = {3e-9, 1e-9, 1e-9};
    const auto likelihoods = echidna::caller::annotate::get_RR_RA_AA_Likelihoods_as_phred_scores(
        std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" ), haplotypeVector, genotypeVector,
        genotypeLikelihoods );

    BOOST_REQUIRE_EQUAL( likelihoods.size(), 3 );  // In ploidy = 2 expect 2 + 1 likelihoods.

    // 1.0, 0.66666, 0.333333 as scaling --> 0, 1.7609, 4.7712
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[0] ), 0 );
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[1] ), 2 );
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[2] ), 5 );
}

BOOST_AUTO_TEST_CASE( testGenotypeLikelihoodsComputationLargeLikelihoods )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );

    varPtr_t snp = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );
    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );
    haplotypeVector.push_back( variantSet_t{} );
    haplotypeVector.push_back( variantSet_t{snp} );

    const auto ploidy = 2;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1}, {} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 3 );
    const std::vector< double > genotypeLikelihoods = {3e9, 1e9, 1e9};
    const auto likelihoods = echidna::caller::annotate::get_RR_RA_AA_Likelihoods_as_phred_scores(
        std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" ), haplotypeVector, genotypeVector,
        genotypeLikelihoods );

    BOOST_REQUIRE_EQUAL( likelihoods.size(), 3 );  // In ploidy = 2 expect 2 + 1 likelihoods.

    // 1.0, 0.66666, 0.333333 as scaling --> 0, 1.7609, 4.7712
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[0] ), 0 );
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[1] ), 2 );
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[2] ), 5 );
}

BOOST_AUTO_TEST_CASE( testGenotypeLikelihoodsComputationForHomAltWithMultipleVariants )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 2 ), "AA" );

    const auto snp0 = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );
    const auto snp1 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "C" );

    HaplotypeVector haplotypeVector( Region( "1", 0, 2 ), referenceSequence );
    auto regionReference = haplotypeVector.regionReferenceSequence();
    haplotypeVector.push_back( variantSet_t{} );            //  hap0
    haplotypeVector.push_back( variantSet_t{snp0} );        //  hap1
    haplotypeVector.push_back( variantSet_t{snp0, snp1} );  //  hap2

    const auto ploidy = 2;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1, 2}, {} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 6 );

    const std::vector< double > genotypeLikelihoods = {1, 1, 1, 1, 1, 1};
    const auto likelihoods = echidna::caller::annotate::get_RR_RA_AA_Likelihoods_as_phred_scores(
        snp0, haplotypeVector, genotypeVector, genotypeLikelihoods );

    BOOST_REQUIRE_EQUAL( likelihoods.size(), 3 );  // In ploidy = 2 expect 2 + 1 likelihoods.

    // Expect 1, 1 * 2 + 1 * 2, 1 + 2 * 1 + 1
    // This is equivalent to 1/4, 1, 1 --> 6.02, 0, 0
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[0] ), 6 );
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[1] ), 0 );
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[2] ), 0 );
}

BOOST_AUTO_TEST_CASE( testGenotypeLikelihoodsComputationForHetCallWithMultipleVariants )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 2 ), "AA" );

    const auto snp0 = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );
    const auto snp1 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "C" );

    HaplotypeVector haplotypeVector( Region( "1", 0, 2 ), referenceSequence );
    haplotypeVector.push_back( variantSet_t{} );            //  hap0
    haplotypeVector.push_back( variantSet_t{snp0} );        //  hap1
    haplotypeVector.push_back( variantSet_t{snp0, snp1} );  //  hap2

    const auto ploidy = 2;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1, 2}, {} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 6 );

    const std::vector< double > genotypeLikelihoods = {1, 1, 1, 1, 1, 1};
    const auto likelihoods = echidna::caller::annotate::get_RR_RA_AA_Likelihoods_as_phred_scores(
        snp1, haplotypeVector, genotypeVector, genotypeLikelihoods );

    BOOST_REQUIRE_EQUAL( likelihoods.size(), 3 );  // In ploidy = 2 expect 2 + 1 likelihoods.

    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[0] ), 0 );
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[1] ), 0 );
    BOOST_CHECK_EQUAL( echidna::stats::roundPhred( likelihoods[2] ), 6 );
}
