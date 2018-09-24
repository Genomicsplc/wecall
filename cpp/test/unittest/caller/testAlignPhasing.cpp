#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "caller/annotation.hpp"
#include "caller/diploid/diploidAnnotate.hpp"
#include "caller/alignPhasing.hpp"
#include "variant/type/variant.hpp"
#include "caller/callSet.hpp"
#include "io/read.hpp"
#include "io/readIntervalTree.hpp"
#include "io/readRange.hpp"

using echidna::io::Read;
using echidna::alignment::Cigar;
using echidna::caller::SetRegions;
using echidna::caller::Region;
using echidna::caller::Call;
using echidna::caller::callVector_t;
using echidna::io::RegionsReads;
using echidna::io::ReadDataset;
using echidna::utils::BasePairSequence;
using echidna::variant::Variant;
using echidna::variant::variantSet_t;
using echidna::caller::Annotation;

// ---------------------------------------------------------------------------------------------------------------
// test: generateVariantSetsFromCalls

BOOST_AUTO_TEST_CASE( shouldGenerateVariantSetFromCallsForEmptyCase )
{
    const auto variantSet1 = echidna::caller::generateVariantSetsFromCalls( {}, 0, 0 );
    BOOST_REQUIRE_EQUAL( variantSet1.size(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldGenerateVariantSetFromCallsForSimpleCase )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'A' ) );
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 2, 3 ), "C" );
    Call call1( variant1, variant1->interval(), 10, 1, {{Call::REF, Call::VAR}} );

    const auto variantSets = echidna::caller::generateVariantSetsFromCalls( {call1}, 0, 2 );

    BOOST_REQUIRE_EQUAL( variantSets.size(), 2 );
    BOOST_REQUIRE_EQUAL( variantSets[0].size(), 0 );
    BOOST_REQUIRE_EQUAL( variantSets[1].size(), 1 );
    BOOST_CHECK_EQUAL( *( variantSets[1].cbegin() ), variant1 );
}

BOOST_AUTO_TEST_CASE( shouldGenerateVariantSetFromCallsForOneSample )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'A' ) );
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 2, 3 ), "C" );
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 3, 4 ), "T" );
    echidna::variant::varPtr_t variant3 = std::make_shared< Variant >( refSequence, Region( "1", 4, 5 ), "G" );
    Call call1( variant1, variant1->interval(), 10, 1, {{Call::REF, Call::VAR}} );
    Call call2( variant2, variant2->interval(), 10, 1, {{Call::VAR, Call::VAR}} );
    Call call3( variant3, variant3->interval(), 10, 1, {{Call::VAR, Call::REF}} );

    const auto ploidy = 2;
    const auto variantSets = echidna::caller::generateVariantSetsFromCalls( {call1, call2, call3}, 0, ploidy );

    BOOST_REQUIRE_EQUAL( variantSets.size(), ploidy );
    BOOST_REQUIRE_EQUAL( variantSets[0].size(), 2 );
    BOOST_CHECK_EQUAL( *( variantSets[0].cbegin() ), variant2 );
    BOOST_CHECK_EQUAL( *( ++variantSets[0].cbegin() ), variant3 );

    BOOST_REQUIRE_EQUAL( variantSets[1].size(), 2 );
    BOOST_CHECK_EQUAL( *( variantSets[1].cbegin() ), variant1 );
    BOOST_CHECK_EQUAL( *( ++variantSets[1].cbegin() ), variant2 );
}

BOOST_AUTO_TEST_CASE( shouldGenerateVariantSetFromCallsForTwoSamples )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'A' ) );
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 2, 3 ), "C" );
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 3, 4 ), "T" );
    Call call1( variant1, variant1->interval(), 10, 2, {{Call::REF, Call::REF}, {Call::REF, Call::VAR}} );
    Call call2( variant2, variant1->interval(), 10, 2, {{Call::VAR, Call::VAR}, {Call::VAR, Call::REF}} );

    const auto ploidy = 2;
    const auto variantSetsSample1 = echidna::caller::generateVariantSetsFromCalls( {call1, call2}, 0, ploidy );

    BOOST_REQUIRE_EQUAL( variantSetsSample1.size(), ploidy );
    BOOST_REQUIRE_EQUAL( variantSetsSample1[0].size(), 1 );
    BOOST_CHECK_EQUAL( *( variantSetsSample1[0].cbegin() ), variant2 );
    BOOST_REQUIRE_EQUAL( variantSetsSample1[1].size(), 1 );
    BOOST_CHECK_EQUAL( *( variantSetsSample1[1].cbegin() ), variant2 );

    const auto variantSetsSample2 = echidna::caller::generateVariantSetsFromCalls( {call1, call2}, 1, ploidy );

    BOOST_REQUIRE_EQUAL( variantSetsSample2.size(), 2 );
    BOOST_REQUIRE_EQUAL( variantSetsSample2[0].size(), 1 );
    BOOST_CHECK_EQUAL( *( variantSetsSample2[0].cbegin() ), variant2 );
    BOOST_REQUIRE_EQUAL( variantSetsSample2[1].size(), 1 );
    BOOST_CHECK_EQUAL( *( variantSetsSample2[1].cbegin() ), variant1 );
}

BOOST_AUTO_TEST_CASE( shouldGenerateVariantSetFromCallsForHaploid )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'A' ) );
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 2, 3 ), "C" );

    Call call1( variant1, variant1->interval(), 10, 1, {{Call::VAR}} );

    const auto ploidy = 1;
    const auto variantSetsSample1 = echidna::caller::generateVariantSetsFromCalls( {call1}, 0, ploidy );

    BOOST_REQUIRE_EQUAL( variantSetsSample1.size(), ploidy );
    BOOST_REQUIRE_EQUAL( variantSetsSample1[0].size(), 1 );
    BOOST_CHECK_EQUAL( *( variantSetsSample1[0].cbegin() ), variant1 );
}

BOOST_AUTO_TEST_CASE( shouldGenerateVariantSetFromCallsForTriploid )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'A' ) );
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 2, 3 ), "C" );
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 3, 4 ), "T" );
    echidna::variant::varPtr_t variant3 = std::make_shared< Variant >( refSequence, Region( "1", 4, 5 ), "G" );
    Call call1( variant1, variant1->interval(), 10, 1, {{Call::REF, Call::VAR, Call::VAR}} );
    Call call2( variant2, variant2->interval(), 10, 1, {{Call::VAR, Call::VAR, Call::REF}} );
    Call call3( variant3, variant3->interval(), 10, 1, {{Call::VAR, Call::REF, Call::REF}} );

    const auto ploidy = 3;
    const auto variantSets = echidna::caller::generateVariantSetsFromCalls( {call1, call2, call3}, 0, ploidy );

    BOOST_REQUIRE_EQUAL( variantSets.size(), ploidy );
    BOOST_REQUIRE_EQUAL( variantSets[0].size(), 2 );
    BOOST_CHECK_EQUAL( *( variantSets[0].cbegin() ), variant2 );
    BOOST_CHECK_EQUAL( *( ++variantSets[0].cbegin() ), variant3 );

    BOOST_REQUIRE_EQUAL( variantSets[1].size(), 2 );
    BOOST_CHECK_EQUAL( *( variantSets[1].cbegin() ), variant1 );
    BOOST_CHECK_EQUAL( *( ++variantSets[1].cbegin() ), variant2 );

    BOOST_REQUIRE_EQUAL( variantSets[2].size(), 1 );
    BOOST_CHECK_EQUAL( *( variantSets[2].cbegin() ), variant1 );
}

// ---------------------------------------------------------------------------------------------------------------
// test: generateHaplotypesFromPhasedClusters

BOOST_AUTO_TEST_CASE( shouldGenerateHaplotypesForClustersWithHetVariants )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'A' ) );
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 2, 3 ), "C" );
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 3, 4 ), "T" );

    // het variants1 (REF, VAR)
    std::vector< echidna::variant::variantSet_t > variantsCluster1;
    variantsCluster1.push_back( {} );
    variantsCluster1.push_back( {variant1} );

    // het variants2 (VAR, REF)
    std::vector< echidna::variant::variantSet_t > variantsCluster2;
    variantsCluster2.push_back( {variant2} );
    variantsCluster2.push_back( {} );

    auto haplotypes = echidna::caller::generateHaplotypesFromPhasedClusters( variantsCluster1, variantsCluster2, 2,
                                                                             refSequence, refSequence->region() );

    BOOST_REQUIRE_EQUAL( haplotypes.size(), 4 );

    BOOST_REQUIRE_EQUAL( haplotypes[0].getVariants().size(), 0 );
    variantSet_t expVariants = {};
    BOOST_CHECK( haplotypes[0].getVariants() == expVariants );
    BOOST_CHECK_EQUAL( haplotypes[0].getId(), 1 );

    BOOST_REQUIRE_EQUAL( haplotypes[1].getVariants().size(), 1 );
    expVariants = {variant2};
    BOOST_CHECK( haplotypes[1].getVariants() == expVariants );
    BOOST_CHECK_EQUAL( haplotypes[1].getId(), 0 );

    BOOST_REQUIRE_EQUAL( haplotypes[2].getVariants().size(), 1 );
    expVariants = {variant1};
    BOOST_CHECK( haplotypes[2].getVariants() == expVariants );
    BOOST_CHECK_EQUAL( haplotypes[2].getId(), 3 );

    BOOST_REQUIRE_EQUAL( haplotypes[3].getVariants().size(), 2 );
    expVariants = {variant1, variant2};
    BOOST_CHECK( haplotypes[3].getVariants() == expVariants );
    BOOST_CHECK_EQUAL( haplotypes[3].getId(), 2 );
}

BOOST_AUTO_TEST_CASE( shouldGenerateIgnoreInvalidHaplotypes )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );
    echidna::variant::varPtr_t large_deletion = std::make_shared< Variant >( refSequence, Region( "1", 2, 10 ), "" );
    echidna::variant::varPtr_t snp = std::make_shared< Variant >( refSequence, Region( "1", 3, 4 ), "T" );

    // hom large deletion (VAR, VAR)
    std::vector< echidna::variant::variantSet_t > variantsCluster1;
    variantsCluster1.push_back( {large_deletion} );
    variantsCluster1.push_back( {large_deletion} );

    // hom snp (VAR, REF)
    std::vector< echidna::variant::variantSet_t > variantsCluster2;
    variantsCluster2.push_back( {snp} );
    variantsCluster2.push_back( {snp} );

    auto haplotypes = echidna::caller::generateHaplotypesFromPhasedClusters( variantsCluster1, variantsCluster2, 2,
                                                                             refSequence, refSequence->region() );

    BOOST_REQUIRE_EQUAL( haplotypes.size(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldGenerateHaplotypesForClustersWithMixedVariants )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'A' ) );
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 2, 3 ), "C" );
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 3, 4 ), "T" );

    // hom variants1
    std::vector< echidna::variant::variantSet_t > variantsCluster1;
    variantsCluster1.push_back( {variant1} );
    variantsCluster1.push_back( {variant1} );

    // het variant2
    std::vector< echidna::variant::variantSet_t > variantsCluster2;
    variantsCluster2.push_back( {variant2} );
    variantsCluster2.push_back( {} );

    auto haplotypes = echidna::caller::generateHaplotypesFromPhasedClusters( variantsCluster1, variantsCluster2, 2,
                                                                             refSequence, refSequence->region() );

    BOOST_REQUIRE_EQUAL( haplotypes.size(), 2 );

    BOOST_REQUIRE_EQUAL( haplotypes[0].getVariants().size(), 1 );
    variantSet_t expVariants = {variant1};
    BOOST_CHECK( haplotypes[0].getVariants() == expVariants );
    BOOST_CHECK_EQUAL( haplotypes[0].getId(), 1 );

    BOOST_REQUIRE_EQUAL( haplotypes[1].getVariants().size(), 2 );
    expVariants = {variant1, variant2};
    BOOST_CHECK( haplotypes[1].getVariants() == expVariants );
    BOOST_CHECK_EQUAL( haplotypes[1].getId(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldGenerateHaplotypesForClustersWithHomVariants )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'A' ) );
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 2, 3 ), "C" );
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 3, 4 ), "T" );

    // hom variants1
    std::vector< echidna::variant::variantSet_t > variantsCluster1;
    variantsCluster1.push_back( {variant1} );
    variantsCluster1.push_back( {variant1} );

    // hom variant2
    std::vector< echidna::variant::variantSet_t > variantsCluster2;
    variantsCluster2.push_back( {variant2} );
    variantsCluster2.push_back( {variant2} );

    auto haplotypes = echidna::caller::generateHaplotypesFromPhasedClusters( variantsCluster1, variantsCluster2, 2,
                                                                             refSequence, refSequence->region() );

    BOOST_REQUIRE_EQUAL( haplotypes.size(), 1 );

    BOOST_REQUIRE_EQUAL( haplotypes[0].getVariants().size(), 2 );
    variantSet_t expVariants = {variant1, variant2};
    BOOST_CHECK( haplotypes[0].getVariants() == expVariants );
    BOOST_CHECK_EQUAL( haplotypes[0].getId(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldGenerateHaplotypesForClustersWithMultipleVariants )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );

    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 2, 3 ), "C" );
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 3, 4 ), "T" );
    echidna::variant::varPtr_t variant3 = std::make_shared< Variant >( refSequence, Region( "1", 5, 6 ), "G" );
    echidna::variant::varPtr_t variant4 = std::make_shared< Variant >( refSequence, Region( "1", 7, 8 ), "C" );

    Call call1( variant1, variant1->interval(), 100, 1, {{Call::REF, Call::VAR}} );
    Call call2( variant2, variant2->interval(), 100, 1, {{Call::VAR, Call::VAR}} );
    Call call3( variant3, variant3->interval(), 100, 1, {{Call::REF, Call::VAR}} );
    Call call4( variant4, variant4->interval(), 100, 1, {{Call::VAR, Call::REF}} );

    // het variants1, hom variant2
    std::vector< echidna::variant::variantSet_t > variantsCluster1;
    variantsCluster1.push_back( {variant2} );
    variantsCluster1.push_back( {variant1, variant2} );

    // het variant3, het variant4
    std::vector< echidna::variant::variantSet_t > variantsCluster2;
    variantsCluster2.push_back( {variant4} );
    variantsCluster2.push_back( {variant3} );

    auto haplotypes = echidna::caller::generateHaplotypesFromPhasedClusters( variantsCluster1, variantsCluster2, 2,
                                                                             refSequence, refSequence->region() );

    BOOST_REQUIRE_EQUAL( haplotypes.size(), 4 );

    BOOST_REQUIRE_EQUAL( haplotypes[0].getVariants().size(), 2 );
    variantSet_t expVariants = {variant2, variant4};
    BOOST_CHECK( haplotypes[0].getVariants() == expVariants );
    BOOST_CHECK_EQUAL( haplotypes[0].getId(), 0 );

    BOOST_REQUIRE_EQUAL( haplotypes[1].getVariants().size(), 2 );
    expVariants = {variant2, variant3};
    BOOST_CHECK( haplotypes[1].getVariants() == expVariants );
    BOOST_CHECK_EQUAL( haplotypes[1].getId(), 1 );

    BOOST_REQUIRE_EQUAL( haplotypes[2].getVariants().size(), 3 );
    expVariants = {variant1, variant2, variant4};
    BOOST_CHECK( haplotypes[2].getVariants() == expVariants );
    BOOST_CHECK_EQUAL( haplotypes[2].getId(), 2 );

    BOOST_REQUIRE_EQUAL( haplotypes[3].getVariants().size(), 3 );
    expVariants = {variant1, variant2, variant3};
    BOOST_CHECK( haplotypes[3].getVariants() == expVariants );
    BOOST_CHECK_EQUAL( haplotypes[3].getId(), 3 );
}

BOOST_AUTO_TEST_CASE( shouldGenerateHaplotypesForHaploid )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'A' ) );
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 2, 3 ), "C" );
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 3, 4 ), "T" );

    // variants1
    std::vector< echidna::variant::variantSet_t > variantsCluster1;
    variantsCluster1.push_back( {variant1} );

    // variant2
    std::vector< echidna::variant::variantSet_t > variantsCluster2;
    variantsCluster2.push_back( {variant2} );

    auto haplotypes = echidna::caller::generateHaplotypesFromPhasedClusters( variantsCluster1, variantsCluster2, 1,
                                                                             refSequence, refSequence->region() );

    BOOST_REQUIRE_EQUAL( haplotypes.size(), 1 );

    BOOST_REQUIRE_EQUAL( haplotypes[0].getVariants().size(), 2 );
    variantSet_t expVariants = {variant1, variant2};
    BOOST_CHECK( haplotypes[0].getVariants() == expVariants );
    BOOST_CHECK_EQUAL( haplotypes[0].getId(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldGenerateHaplotypesForClustersWithHetVariantsForTriploid )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'A' ) );
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 2, 3 ), "C" );
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 3, 4 ), "T" );

    // 0/0/1 variants1
    std::vector< echidna::variant::variantSet_t > variantsCluster1;
    variantsCluster1.push_back( {} );
    variantsCluster1.push_back( {} );
    variantsCluster1.push_back( {variant1} );

    // 0/1/1 variant2
    std::vector< echidna::variant::variantSet_t > variantsCluster2;
    variantsCluster2.push_back( {} );
    variantsCluster2.push_back( {variant2} );
    variantsCluster2.push_back( {variant2} );

    auto haplotypes = echidna::caller::generateHaplotypesFromPhasedClusters( variantsCluster1, variantsCluster2, 3,
                                                                             refSequence, refSequence->region() );

    BOOST_REQUIRE_EQUAL( haplotypes.size(), 4 );

    BOOST_REQUIRE_EQUAL( haplotypes[0].getVariants().size(), 0 );
    variantSet_t expVariants = {};
    BOOST_CHECK( haplotypes[0].getVariants() == expVariants );
    BOOST_CHECK_EQUAL( haplotypes[0].getId(), 0 );

    BOOST_REQUIRE_EQUAL( haplotypes[1].getVariants().size(), 1 );
    expVariants = {variant2};
    BOOST_CHECK( haplotypes[1].getVariants() == expVariants );
    BOOST_CHECK_EQUAL( haplotypes[1].getId(), 1 );

    BOOST_REQUIRE_EQUAL( haplotypes[2].getVariants().size(), 1 );
    expVariants = {variant1};
    BOOST_CHECK( haplotypes[2].getVariants() == expVariants );
    BOOST_CHECK_EQUAL( haplotypes[2].getId(), 6 );

    BOOST_REQUIRE_EQUAL( haplotypes[3].getVariants().size(), 2 );
    expVariants = {variant1, variant2};
    BOOST_CHECK( haplotypes[3].getVariants() == expVariants );
    BOOST_CHECK_EQUAL( haplotypes[3].getId(), 7 );
}

// ---------------------------------------------------------------------------------------------------------------
// test: alignPhasingOfCallsBasedOnGenotypeAndMergedHaplotypes

BOOST_AUTO_TEST_CASE( shouldAllignCallsForGenotypeWithHetVariantsOnSameRead )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );

    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 2, 3 ), "C" );
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 6, 7 ), "T" );

    const int64_t phaseId = 101;
    Call call2( variant2, variant2->interval(), 100, 1, {{Call::VAR, Call::REF}} );
    call2.samples[0].addAnnotation( Annotation::PS, phaseId + 1 );

    // het variants1 (REF, VAR)
    std::vector< echidna::variant::variantSet_t > variantsCluster1;
    variantsCluster1.push_back( {} );
    variantsCluster1.push_back( {variant1} );

    echidna::variant::HaplotypeVector combinedHaplotypes( refSequence->region(), refSequence );
    combinedHaplotypes.push_back( {}, 2 );
    combinedHaplotypes.push_back( {variant2}, 3 );
    combinedHaplotypes.push_back( {variant1}, 0 );
    combinedHaplotypes.push_back( {variant1, variant2}, 1 );

    // calling 1x haplotype0 aqnd 1x haplotype3
    auto genotype = echidna::variant::Genotype( echidna::variant::haplotypeAndCount_t{{0, 1}, {3, 1}} );
    echidna::variant::genotypePtr_t genotypePtr = std::make_shared< echidna::variant::Genotype >( genotype );

    callVector_t calls = {call2};
    echidna::caller::alignPhasingOfCallsBasedOnGenotypeAndMergedHaplotypes( genotypePtr, 2, 0, phaseId,
                                                                            combinedHaplotypes, false, calls );

    BOOST_REQUIRE_EQUAL( calls.size(), 1 );
    BOOST_REQUIRE_EQUAL( calls[0].samples.size(), 1 );
    BOOST_REQUIRE_EQUAL( calls[0].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( calls[0].samples[0].genotypeCalls[0], Call::REF );
    BOOST_CHECK_EQUAL( calls[0].samples[0].genotypeCalls[1], Call::VAR );
    BOOST_CHECK_EQUAL( calls[0].samples[0].getAnnotation( Annotation::PS ), phaseId );
}

BOOST_AUTO_TEST_CASE( shouldAllignCallsForGenotypeWithHetVariantsOnDifferentRead )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );

    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 2, 3 ), "C" );
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 6, 7 ), "T" );

    const int64_t phaseId = 101;
    Call call2( variant2, variant2->interval(), 100, 1, {{Call::REF, Call::VAR}} );
    call2.samples[0].addAnnotation( Annotation::PS, phaseId + 1 );

    // het variants1 (REF, VAR)
    std::vector< echidna::variant::variantSet_t > variantsCluster1;
    variantsCluster1.push_back( {} );
    variantsCluster1.push_back( {variant1} );

    echidna::variant::HaplotypeVector combinedHaplotypes( refSequence->region(), refSequence );
    combinedHaplotypes.push_back( {}, 0 );
    combinedHaplotypes.push_back( {variant2}, 1 );
    combinedHaplotypes.push_back( {variant1}, 2 );
    combinedHaplotypes.push_back( {variant1, variant2}, 3 );

    // calling 1x haplotype1 aqnd 1x haplotype2
    auto genotype = echidna::variant::Genotype( echidna::variant::haplotypeAndCount_t{{1, 1}, {2, 1}} );
    echidna::variant::genotypePtr_t genotypePtr = std::make_shared< echidna::variant::Genotype >( genotype );

    callVector_t calls = {call2};
    echidna::caller::alignPhasingOfCallsBasedOnGenotypeAndMergedHaplotypes( genotypePtr, 2, 0, phaseId,
                                                                            combinedHaplotypes, false, calls );

    BOOST_REQUIRE_EQUAL( calls.size(), 1 );
    BOOST_REQUIRE_EQUAL( calls[0].samples.size(), 1 );
    BOOST_REQUIRE_EQUAL( calls[0].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( calls[0].samples[0].genotypeCalls[0], Call::VAR );
    BOOST_CHECK_EQUAL( calls[0].samples[0].genotypeCalls[1], Call::REF );
    BOOST_CHECK_EQUAL( calls[0].samples[0].getAnnotation( Annotation::PS ), phaseId );
}

BOOST_AUTO_TEST_CASE( shouldAllignCallsForGenotypeWithHomVariant )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );

    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 2, 3 ), "C" );
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 6, 7 ), "T" );

    const int64_t phaseId = 101;
    Call call2( variant2, variant2->interval(), 100, 1, {{Call::VAR, Call::VAR}} );
    call2.samples[0].addAnnotation( Annotation::PS, phaseId + 1 );

    // het variants1 (REF, VAR)
    std::vector< echidna::variant::variantSet_t > variantsCluster1;
    variantsCluster1.push_back( {} );
    variantsCluster1.push_back( {variant1} );

    echidna::variant::HaplotypeVector combinedHaplotypes( refSequence->region(), refSequence );
    combinedHaplotypes.push_back( {variant2}, 0 );
    combinedHaplotypes.push_back( {variant1, variant2}, 2 );

    // calling 1x haplotype0 aqnd 1x haplotype1
    auto genotype = echidna::variant::Genotype( echidna::variant::haplotypeAndCount_t{{0, 1}, {1, 1}} );
    echidna::variant::genotypePtr_t genotypePtr = std::make_shared< echidna::variant::Genotype >( genotype );

    callVector_t calls = {call2};
    echidna::caller::alignPhasingOfCallsBasedOnGenotypeAndMergedHaplotypes( genotypePtr, 2, 0, phaseId,
                                                                            combinedHaplotypes, true, calls );

    BOOST_REQUIRE_EQUAL( calls.size(), 1 );
    BOOST_REQUIRE_EQUAL( calls[0].samples.size(), 1 );
    BOOST_REQUIRE_EQUAL( calls[0].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( calls[0].samples[0].genotypeCalls[0], Call::VAR );
    BOOST_CHECK_EQUAL( calls[0].samples[0].genotypeCalls[1], Call::VAR );
    BOOST_CHECK_EQUAL( calls[0].samples[0].getAnnotation( Annotation::PS ), phaseId );
}

BOOST_AUTO_TEST_CASE( shouldAllignCallsForGenotypeForVariantsOnSameReadForSecondSample )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );

    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 2, 3 ), "C" );
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 6, 7 ), "T" );

    const int64_t phaseId = 101;
    Call call2( variant2, variant2->interval(), 100, 2, {{Call::VAR, Call::VAR}, {Call::VAR, Call::REF}} );
    const size_t sampleIndex = 1;
    call2.samples[sampleIndex].addAnnotation( Annotation::PS, phaseId + 1 );

    // het variants1 (REF, VAR) for sample 2
    std::vector< echidna::variant::variantSet_t > variantsCluster1;
    variantsCluster1.push_back( {} );
    variantsCluster1.push_back( {variant1} );

    echidna::variant::HaplotypeVector combinedHaplotypes( refSequence->region(), refSequence );
    combinedHaplotypes.push_back( {}, 2 );
    combinedHaplotypes.push_back( {variant2}, 3 );
    combinedHaplotypes.push_back( {variant1}, 0 );
    combinedHaplotypes.push_back( {variant1, variant2}, 1 );

    // calling 1x haplotype0 aqnd 1x haplotype3
    const auto genotype = echidna::variant::Genotype( echidna::variant::haplotypeAndCount_t{{0, 1}, {3, 1}} );
    echidna::variant::genotypePtr_t genotypePtr = std::make_shared< echidna::variant::Genotype >( genotype );

    callVector_t calls = {call2};
    echidna::caller::alignPhasingOfCallsBasedOnGenotypeAndMergedHaplotypes( genotypePtr, 2, sampleIndex, phaseId,
                                                                            combinedHaplotypes, false, calls );

    BOOST_REQUIRE_EQUAL( calls.size(), 1 );
    BOOST_REQUIRE_EQUAL( calls[0].samples.size(), 2 );
    BOOST_REQUIRE_EQUAL( calls[0].samples[sampleIndex].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( calls[0].samples[sampleIndex].genotypeCalls[0], Call::REF );
    BOOST_CHECK_EQUAL( calls[0].samples[sampleIndex].genotypeCalls[1], Call::VAR );
    BOOST_CHECK_EQUAL( calls[0].samples[sampleIndex].getAnnotation( Annotation::PS ), phaseId );
}

// ---------------------------------------------------------------------------------------------------------------
// test: alignPhasingForSample

BOOST_AUTO_TEST_CASE( shouldAlignPhaseSetsOfClustersForVariantOnDifferentRead )
{
    // two variants  on different strands
    auto reference = "GAGGGTCCTGCAAGGAACTGCGGGAAGTCT";
    // var1          "          .C......            " on strand1
    // var2          "          ...T....            " on strand2

    // generate overlapping reads
    // reads need to be at least kmerLength and reference sequence needs to be sufficiently long for padding
    Region region = Region( "1", 0, 30 );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( region, reference );

    const int64_t startPos = 10;
    const uint8_t mapQual = 10;
    const auto read1 =
        std::make_shared< Read >( BasePairSequence( "CCAGGAAC" ), std::string( 8, 'Q' ), "0", Cigar( "8M" ), 0,
                                  startPos, 0, mapQual, 0, 0, 0, refSequence, "read1" );
    const auto read2 =
        std::make_shared< Read >( BasePairSequence( "CAATGAAC" ), std::string( 8, 'Q' ), "0", Cigar( "8M" ), 0,
                                  startPos, 0, mapQual, 0, 0, 0, refSequence, "read2" );

    echidna::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );
    readContainer.insert( read2 );

    RegionsReads regionSetReads( region, readContainer.getFullRange(), 0 );

    auto sampleName = "sample1";
    echidna::io::perSampleRegionsReads_t overlappingReads = {{sampleName, regionSetReads}};

    // generate call vector for previous cluster
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 11, 12 ), "C" );
    Call call1( variant1, variant1->interval(), 100, 1, {{Call::REF, Call::VAR}} );
    const int64_t phaseId = 101;
    call1.samples[0].addAnnotation( Annotation::PS, phaseId );

    // generate call vector for current cluster
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 13, 14 ), "T" );
    Call call2( variant2, variant2->interval(), 100, 1, {{Call::REF, Call::VAR}} );
    call2.samples[0].addAnnotation( Annotation::PS, phaseId + 1 );

    // align phase sets
    callVector_t alignedCalls = {call2};
    echidna::caller::alignPhasingForSample( alignedCalls, {call1}, 2, 0, sampleName, refSequence, overlappingReads,
                                            Region( "1", 10, 12 ), refSequence->region() );

    BOOST_REQUIRE_EQUAL( alignedCalls.size(), 1 );
    BOOST_REQUIRE_EQUAL( alignedCalls[0].samples.size(), 1 );
    BOOST_REQUIRE_EQUAL( alignedCalls[0].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].genotypeCalls[0], Call::VAR );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].genotypeCalls[1], Call::REF );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].getAnnotation( Annotation::PS ), phaseId );
}

BOOST_AUTO_TEST_CASE( shouldAllignPhaseSetsOfClustersForVariantOnSameRead )
{
    // two variants on same strand
    auto reference = "GAGGGTCCTGCAAGGAACTGCGGGAAGTCT";
    // var1, var2    "          .C.T....            " on strand1
    // no variants   "          ........            " on strand2

    // generate overlapping reads
    // reads need to be at least kmerLength and reference sequence needs to be sufficiently long for padding
    Region region = Region( "1", 0, 30 );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( region, reference );

    const int64_t startPos = 10;
    const uint8_t mapQual = 10;
    const auto read1 =
        std::make_shared< Read >( BasePairSequence( "CCATGAAC" ), std::string( 8, 'Q' ), "0", Cigar( "8M" ), 0,
                                  startPos, 0, mapQual, 0, 0, 0, refSequence, "read1" );
    const auto read2 =
        std::make_shared< Read >( BasePairSequence( "CAAGGAAC" ), std::string( 8, 'Q' ), "0", Cigar( "8M" ), 0,
                                  startPos, 0, mapQual, 0, 0, 0, refSequence, "read2" );

    echidna::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );
    readContainer.insert( read2 );

    RegionsReads regionSetReads( region, readContainer.getFullRange(), 0 );

    auto sampleName = "sample1";
    echidna::io::perSampleRegionsReads_t overlappingReads = {{sampleName, regionSetReads}};

    // generate call vector for previous cluster
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 11, 12 ), "C" );
    Call call1( variant1, variant1->interval(), 100, 1, {{Call::REF, Call::VAR}} );
    const int64_t phaseId = 101;
    call1.samples[0].addAnnotation( Annotation::PS, phaseId );

    // generate call vector for current cluster
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 13, 14 ), "T" );
    Call call2( variant2, variant2->interval(), 100, 1, {{Call::VAR, Call::REF}} );
    call2.samples[0].addAnnotation( Annotation::PS, phaseId + 1 );

    // align phase sets
    callVector_t alignedCalls = {call2};
    echidna::caller::alignPhasingForSample( alignedCalls, {call1}, 2, 0, sampleName, refSequence, overlappingReads,
                                            Region( "1", 10, 12 ), refSequence->region() );

    BOOST_REQUIRE_EQUAL( alignedCalls.size(), 1 );
    BOOST_REQUIRE_EQUAL( alignedCalls[0].samples.size(), 1 );
    BOOST_REQUIRE_EQUAL( alignedCalls[0].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].genotypeCalls[0], Call::REF );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].genotypeCalls[1], Call::VAR );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].getAnnotation( Annotation::PS ), phaseId );
}

BOOST_AUTO_TEST_CASE( shouldAllignPhaseSetsOfClustersForMultipleVariants )
{
    // four variants on different strands
    auto reference = "GAGGGTCCTGCAAGGAACTGCGGGAAGTCT";
    // var2, var4    "       ......T...A...         " on strand 1
    // var1, var3    "       ....C...C.....         " on strand 2

    // generate overlapping reads
    // reads need to be at least kmerLength and reference sequence needs to be sufficiently long for padding
    Region region = Region( "1", 0, 30 );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( region, reference );

    const int64_t startPos = 10;
    const uint8_t mapQual = 10;
    const auto read1 =
        std::make_shared< Read >( BasePairSequence( "CTGCAATGAAATGC" ), std::string( 14, 'Q' ), "0", Cigar( "14M" ), 0,
                                  startPos, 0, mapQual, 0, 0, 0, refSequence, "read1" );
    const auto read2 =
        std::make_shared< Read >( BasePairSequence( "CTGCCATGCACTGC" ), std::string( 14, 'Q' ), "0", Cigar( "14M" ), 0,
                                  startPos, 0, mapQual, 0, 0, 0, refSequence, "read2" );

    echidna::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );
    readContainer.insert( read2 );

    RegionsReads regionSetReads( region, readContainer.getFullRange(), 0 );

    auto sampleName = "sample1";
    echidna::io::perSampleRegionsReads_t overlappingReads = {{sampleName, regionSetReads}};

    // generate call vector for previous cluster
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 11, 12 ), "C" );
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 13, 14 ), "T" );
    Call call1( variant1, variant1->interval(), 100, 1, {{Call::REF, Call::VAR}} );
    Call call2( variant2, variant2->interval(), 100, 1, {{Call::VAR, Call::VAR}} );
    const int64_t phaseId = 101;
    call1.samples[0].addAnnotation( Annotation::PS, phaseId );
    call2.samples[0].addAnnotation( Annotation::PS, phaseId );

    // generate call vector for current cluster
    echidna::variant::varPtr_t variant3 = std::make_shared< Variant >( refSequence, Region( "1", 15, 16 ), "C" );
    echidna::variant::varPtr_t variant4 = std::make_shared< Variant >( refSequence, Region( "1", 17, 18 ), "A" );
    Call call3( variant3, variant3->interval(), 100, 1, {{Call::REF, Call::VAR}} );
    Call call4( variant4, variant4->interval(), 100, 1, {{Call::VAR, Call::REF}} );
    call3.samples[0].addAnnotation( Annotation::PS, phaseId + 1 );
    call4.samples[0].addAnnotation( Annotation::PS, phaseId + 1 );

    // align phase sets
    callVector_t alignedCalls = {call3, call4};
    echidna::caller::alignPhasingForSample( alignedCalls, {call1, call2}, 2, 0, sampleName, refSequence,
                                            overlappingReads, Region( "1", 10, 14 ), refSequence->region() );

    BOOST_REQUIRE_EQUAL( alignedCalls.size(), 2 );
    BOOST_REQUIRE_EQUAL( alignedCalls[0].samples.size(), 1 );
    BOOST_REQUIRE_EQUAL( alignedCalls[0].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].genotypeCalls[0], Call::REF );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].genotypeCalls[1], Call::VAR );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].getAnnotation( Annotation::PS ), phaseId );

    BOOST_REQUIRE_EQUAL( alignedCalls[1].samples.size(), 1 );
    BOOST_REQUIRE_EQUAL( alignedCalls[1].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( alignedCalls[1].samples[0].genotypeCalls[0], Call::VAR );
    BOOST_CHECK_EQUAL( alignedCalls[1].samples[0].genotypeCalls[1], Call::REF );
    BOOST_CHECK_EQUAL( alignedCalls[1].samples[0].getAnnotation( Annotation::PS ), phaseId );
}

//// ---------------------------------------------------------------------------------------------------------------
// test: alignPhasingBetweenClusters

BOOST_AUTO_TEST_CASE( shouldAllignPhaseSetsOfClustersForOneSample )
{
    // two variants  on different strands
    auto reference = "GAGGGTCCTGCAAGGAACTGCGGGAAGTCT";
    // var1          "          .C......            " on strand1
    // var2          "          ...T....            " on strand2

    // generate overlapping reads
    // reads need to be at least kmerLength and reference sequence needs to be sufficiently long for padding
    Region region = Region( "1", 0, 30 );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( region, reference );

    const int64_t startPos = 10;
    const uint8_t mapQual = 10;
    const auto read1 =
        std::make_shared< Read >( BasePairSequence( "CCAGGAAC" ), std::string( 8, 'Q' ), "0", Cigar( "8M" ), 0,
                                  startPos, 0, mapQual, 0, 0, 0, refSequence, "read1" );
    const auto read2 =
        std::make_shared< Read >( BasePairSequence( "CAATGAAC" ), std::string( 8, 'Q' ), "0", Cigar( "8M" ), 0,
                                  startPos, 0, mapQual, 0, 0, 0, refSequence, "read2" );

    echidna::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );
    readContainer.insert( read2 );

    RegionsReads regionSetReads( region, readContainer.getFullRange(), 0 );

    std::vector< std::string > samples = {"sample1"};
    echidna::io::perSampleRegionsReads_t overlappingReads = {{samples.front(), regionSetReads}};

    // generate call vector for previous cluster
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 11, 12 ), "C" );
    Call call1( variant1, variant1->interval(), 100, 1, {{Call::REF, Call::VAR}} );
    const int64_t phaseId = 101;
    call1.samples[0].addAnnotation( Annotation::PS, phaseId );

    // generate call vector for current cluster
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 13, 14 ), "T" );
    Call call2( variant2, variant2->interval(), 100, 1, {{Call::REF, Call::VAR}} );
    call2.samples[0].addAnnotation( Annotation::PS, phaseId + 1 );

    // align phase sets
    callVector_t alignedCalls = {call2};
    echidna::caller::alignPhasingBetweenClusters( alignedCalls, {call1}, Region( "1", 13, 15 ), Region( "1", 10, 12 ),
                                                  region, refSequence, {2}, overlappingReads, samples );

    BOOST_REQUIRE_EQUAL( alignedCalls.size(), 1 );
    BOOST_REQUIRE_EQUAL( alignedCalls[0].samples.size(), 1 );
    BOOST_REQUIRE_EQUAL( alignedCalls[0].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].genotypeCalls[0], Call::VAR );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].genotypeCalls[1], Call::REF );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].getAnnotation( Annotation::PS ), phaseId );
}

BOOST_AUTO_TEST_CASE( shouldNotAllignPhaseSetsOfClustersForOneSampleWithoutOverlappingReads )
{
    // two variants  on different strands
    auto reference = "GAGGGTCCTGCAAGGAACTGCGGGAAGTCT";
    // var1          "  .C......                    " on strand1
    // var2          "                   ...T....   " on strand2

    // generate overlapping reads
    // reads need to be at least kmerLength and reference sequence needs to be sufficiently long for padding
    Region region = Region( "1", 0, 30 );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( region, reference );

    const uint8_t mapQual = 10;
    auto read1 = std::make_shared< Read >( BasePairSequence( "GCGTCCTG" ), std::string( 8, 'Q' ), "0", Cigar( "8M" ), 0,
                                           2, 0, mapQual, 0, 0, 0, refSequence, "read1" );
    const auto read2 = std::make_shared< Read >( BasePairSequence( "GCGTGAAG" ), std::string( 8, 'Q' ), "0",
                                                 Cigar( "8M" ), 0, 19, 0, mapQual, 0, 0, 0, refSequence, "read2" );

    echidna::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );
    readContainer.insert( read2 );

    RegionsReads regionSetReads( region, readContainer.getFullRange(), 0 );

    std::vector< std::string > samples = {"sample1"};
    echidna::io::perSampleRegionsReads_t overlappingReads = {{samples.front(), regionSetReads}};

    // generate call vector for previous cluster
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 2, 3 ), "C" );
    Call call1( variant1, variant1->interval(), 100, 1, {{Call::REF, Call::VAR}} );
    const int64_t phaseId = 101;
    call1.samples[0].addAnnotation( Annotation::PS, phaseId );

    // generate call vector for current cluster
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 19, 20 ), "T" );
    Call call2( variant2, variant2->interval(), 100, 1, {{Call::REF, Call::VAR}} );
    call2.samples[0].addAnnotation( Annotation::PS, phaseId + 1 );

    // align phase sets
    callVector_t alignedCalls = {call2};
    echidna::caller::alignPhasingBetweenClusters( alignedCalls, {call1}, Region( "1", 14, 25 ), Region( "1", 0, 10 ),
                                                  region, refSequence, {2}, overlappingReads, samples );

    BOOST_REQUIRE_EQUAL( alignedCalls.size(), 1 );
    BOOST_REQUIRE_EQUAL( alignedCalls[0].samples.size(), 1 );
    BOOST_REQUIRE_EQUAL( alignedCalls[0].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].genotypeCalls[0], Call::REF );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].genotypeCalls[1], Call::VAR );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].getAnnotation( Annotation::PS ), phaseId + 1 );
}

BOOST_AUTO_TEST_CASE( shouldAllignPhaseSetsOfClustersForOneSampleWithOverlappingVariants )
{
    // two variants  on different strands
    auto reference = "GAGGGTCCTGCAAGGAACTGCGGGAAGTCT";
    // var1, var2    "          .C.T....            " on strand1
    // var3          "          ...A....            " on strand2

    // generate overlapping reads
    // reads need to be at least kmerLength and reference sequence needs to be sufficiently long for padding
    Region region = Region( "1", 0, 30 );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( region, reference );

    const int64_t startPos = 10;
    const uint8_t mapQual = 10;
    const auto read1 =
        std::make_shared< Read >( BasePairSequence( "CCATGAAC" ), std::string( 8, 'Q' ), "0", Cigar( "8M" ), 0,
                                  startPos, 0, mapQual, 0, 0, 0, refSequence, "read1" );

    const auto read2 =
        std::make_shared< Read >( BasePairSequence( "CAAAGAAC" ), std::string( 8, 'Q' ), "0", Cigar( "8M" ), 0,
                                  startPos, 0, mapQual, 0, 0, 0, refSequence, "read2" );

    echidna::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );
    readContainer.insert( read2 );

    RegionsReads regionSetReads( region, readContainer.getFullRange(), 0 );

    std::vector< std::string > samples = {"sample1"};
    echidna::io::perSampleRegionsReads_t overlappingReads = {{samples.front(), regionSetReads}};

    // generate call vector for previous cluster
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 11, 12 ), "C" );
    Call call1( variant1, variant1->interval(), 100, 1, {{Call::REF, Call::VAR}} );
    const int64_t phaseId = 101;
    call1.samples[0].addAnnotation( Annotation::PS, phaseId );

    // generate call vector for current cluster
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 13, 14 ), "T" );
    echidna::variant::varPtr_t variant3 = std::make_shared< Variant >( refSequence, Region( "1", 13, 14 ), "A" );
    Call call2( variant2, variant2->interval(), 100, 1, {{Call::VAR, Call::UNKNOWN}} );
    Call call3( variant3, variant3->interval(), 100, 1, {{Call::UNKNOWN, Call::VAR}} );
    call2.samples[0].addAnnotation( Annotation::PS, phaseId + 1 );
    call3.samples[0].addAnnotation( Annotation::PS, phaseId + 1 );

    // align phase sets
    callVector_t alignedCalls = {call2, call3};
    echidna::caller::alignPhasingBetweenClusters( alignedCalls, {call1}, Region( "1", 13, 15 ), Region( "1", 10, 12 ),
                                                  region, refSequence, {2}, overlappingReads, samples );

    BOOST_REQUIRE_EQUAL( alignedCalls.size(), 2 );
    BOOST_REQUIRE_EQUAL( alignedCalls[0].samples.size(), 1 );
    BOOST_REQUIRE_EQUAL( alignedCalls[0].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].genotypeCalls[0], Call::UNKNOWN );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].genotypeCalls[1], Call::VAR );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].getAnnotation( Annotation::PS ), phaseId );

    BOOST_REQUIRE_EQUAL( alignedCalls[1].samples.size(), 1 );
    BOOST_REQUIRE_EQUAL( alignedCalls[1].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( alignedCalls[1].samples[0].genotypeCalls[0], Call::VAR );
    BOOST_CHECK_EQUAL( alignedCalls[1].samples[0].genotypeCalls[1], Call::UNKNOWN );
    BOOST_CHECK_EQUAL( alignedCalls[1].samples[0].getAnnotation( Annotation::PS ), phaseId );
}

BOOST_AUTO_TEST_CASE( shouldAllignPhaseSetsOfClustersForTwoSamplesWithOverlappingVariants )
{
    // two variants  on different strands
    auto reference = "GAGGGTCCTGCAAGGAACTGCGGGAAGTCT";
    // var1, var2    "          .C.A....            " on strand1 of sample 1
    // var2          "          ...A....            " on strand2 of sample 1
    // var1, var2    "          .C.A....            " on strand1 of sample 1
    // var3          "          ...T....            " on strand2 of sample 1

    // generate overlapping reads
    // reads need to be at least kmerLength and reference sequence needs to be sufficiently long for padding
    Region region = Region( "1", 0, 30 );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( region, reference );

    const int64_t startPos = 10;
    const uint8_t mapQual = 10;
    const auto read1 =
        std::make_shared< Read >( BasePairSequence( "CCAAGAAC" ), std::string( 8, 'Q' ), "0", Cigar( "8M" ), 0,
                                  startPos, 0, mapQual, 0, 0, 0, refSequence, "read1" );

    const auto read2 =
        std::make_shared< Read >( BasePairSequence( "CAAAGAAC" ), std::string( 8, 'Q' ), "0", Cigar( "8M" ), 0,
                                  startPos, 0, mapQual, 0, 0, 0, refSequence, "read2" );

    const auto read3 =
        std::make_shared< Read >( BasePairSequence( "CCAAGAAC" ), std::string( 8, 'Q' ), "0", Cigar( "8M" ), 0,
                                  startPos, 0, mapQual, 0, 0, 0, refSequence, "read3" );

    const auto read4 =
        std::make_shared< Read >( BasePairSequence( "CAATGAAC" ), std::string( 8, 'Q' ), "0", Cigar( "8M" ), 0,
                                  startPos, 0, mapQual, 0, 0, 0, refSequence, "read4" );

    echidna::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );
    readContainer.insert( read2 );
    readContainer.insert( read3 );
    readContainer.insert( read4 );

    const auto nSamples = 2;

    RegionsReads regionSetReads( region, readContainer.getFullRange(), 0 );

    std::vector< std::string > samples = {"sample1", "sample2"};
    echidna::io::perSampleRegionsReads_t overlappingReads = {{samples[0], regionSetReads},
                                                             {samples[1], regionSetReads}};

    // generate call vector for previous cluster
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, Region( "1", 11, 12 ), "C" );
    Call call1( variant1, variant1->interval(), 100, nSamples, {{Call::REF, Call::VAR}, {Call::REF, Call::VAR}} );
    const int64_t phaseId = 101;
    call1.samples[0].addAnnotation( Annotation::PS, phaseId );
    call1.samples[1].addAnnotation( Annotation::PS, phaseId );

    // generate call vector for current cluster
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequence, Region( "1", 13, 14 ), "A" );
    echidna::variant::varPtr_t variant3 = std::make_shared< Variant >( refSequence, Region( "1", 13, 14 ), "T" );
    Call call2( variant2, variant2->interval(), 100, nSamples, {{Call::VAR, Call::VAR}, {Call::VAR, Call::UNKNOWN}} );
    Call call3( variant3, variant3->interval(), 100, nSamples,
                {{Call::UNKNOWN, Call::UNKNOWN}, {Call::UNKNOWN, Call::VAR}} );
    call2.samples[0].addAnnotation( Annotation::PS, phaseId + 1 );
    call2.samples[1].addAnnotation( Annotation::PS, phaseId + 1 );
    call3.samples[0].addAnnotation( Annotation::PS, phaseId + 1 );
    call3.samples[1].addAnnotation( Annotation::PS, phaseId + 1 );

    // align phase sets
    callVector_t alignedCalls = {call2, call3};
    echidna::caller::alignPhasingBetweenClusters( alignedCalls, {call1}, Region( "1", 13, 15 ), Region( "1", 10, 12 ),
                                                  region, refSequence, {2, 2}, overlappingReads, samples );
    // var2 calls (1/1 for sample1, ./1 for sample2)
    BOOST_REQUIRE_EQUAL( alignedCalls.size(), 2 );
    BOOST_REQUIRE_EQUAL( alignedCalls[0].samples.size(), 2 );
    BOOST_REQUIRE_EQUAL( alignedCalls[0].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].genotypeCalls[0], Call::VAR );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].genotypeCalls[1], Call::VAR );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[0].getAnnotation( Annotation::PS ), phaseId );

    BOOST_REQUIRE_EQUAL( alignedCalls[0].samples[1].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[1].genotypeCalls[0], Call::UNKNOWN );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[1].genotypeCalls[1], Call::VAR );
    BOOST_CHECK_EQUAL( alignedCalls[0].samples[1].getAnnotation( Annotation::PS ), phaseId );

    // var3 calls (0/0 for sample1, 1/. for sample2)
    BOOST_REQUIRE_EQUAL( alignedCalls[1].samples.size(), 2 );
    BOOST_REQUIRE_EQUAL( alignedCalls[1].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( alignedCalls[1].samples[0].genotypeCalls[0], Call::UNKNOWN );
    BOOST_CHECK_EQUAL( alignedCalls[1].samples[0].genotypeCalls[1], Call::UNKNOWN );
    BOOST_CHECK_EQUAL( alignedCalls[1].samples[0].getAnnotation( Annotation::PS ), phaseId );

    BOOST_REQUIRE_EQUAL( alignedCalls[1].samples[1].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( alignedCalls[1].samples[1].genotypeCalls[0], Call::VAR );
    BOOST_CHECK_EQUAL( alignedCalls[1].samples[1].genotypeCalls[1], Call::UNKNOWN );
    BOOST_CHECK_EQUAL( alignedCalls[1].samples[1].getAnnotation( Annotation::PS ), phaseId );
}
