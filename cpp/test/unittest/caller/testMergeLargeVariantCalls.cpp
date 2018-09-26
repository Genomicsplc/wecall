#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "variant/type/variant.hpp"
#include "caller/callSet.hpp"
#include "caller/mergeLargeVariantCalls.hpp"

using echidna::caller::Call;
using echidna::variant::Variant;
using echidna::caller::Region;
using echidna::caller::CallComp;
using echidna::caller::callVector_t;
using echidna::utils::Interval;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( shouldFilterReferenceCallsForOneSample )
{
    auto refSequenceChr =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'A' ) );

    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequenceChr, Region( "1", 1, 2 ), "C" );

    Call call( variant1, variant1->interval(), 10, 1, {{Call::REF, Call::VAR}} );
    Call refCall( nullptr, Interval( 2, 5 ), 10, 1, {{}} );

    const auto callMerger = echidna::caller::lvcMerger();

    // remove ref call when ploidy < defaultPloidy
    auto result = callMerger.removeReferenceCalls( {call, refCall}, {1}, {2} );
    BOOST_CHECK_EQUAL( result.size(), 1 );
    BOOST_CHECK_EQUAL( result[0].var, variant1 );

    // remove ref call when ploidy == defaultPloidy
    result = callMerger.removeReferenceCalls( {call, refCall}, {2}, {2} );
    BOOST_CHECK_EQUAL( result.size(), 2 );
    BOOST_CHECK_EQUAL( result[0].var, variant1 );
    BOOST_CHECK( result[1].isRefCall() );
}

BOOST_AUTO_TEST_CASE( shouldFilterReferenceCallsForTwoSamples )
{
    auto refSequenceChr =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'A' ) );

    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequenceChr, Region( "1", 1, 2 ), "C" );

    Call call( variant1, variant1->interval(), 10, 2, {{Call::REF, Call::VAR}, {Call::REF, Call::REF}} );
    Call refCall( nullptr, Interval( 2, 5 ), 10, 2, {{}, {}} );

    const auto callMerger = echidna::caller::lvcMerger();

    // remove ref call when ploidy < defaultPloidy for at least one sample
    auto result = callMerger.removeReferenceCalls( {call, refCall}, {1, 2}, {2, 2} );
    BOOST_CHECK_EQUAL( result.size(), 1 );
    BOOST_CHECK_EQUAL( result[0].var, variant1 );

    // remove ref call when ploidy == defaultPloidy for all samples
    result = callMerger.removeReferenceCalls( {call, refCall}, {2, 2}, {2, 2} );
    BOOST_CHECK_EQUAL( result.size(), 2 );
    BOOST_CHECK_EQUAL( result[0].var, variant1 );
    BOOST_CHECK( result[1].isRefCall() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( shouldRetrieveOverlappingVariant )
{
    auto region1 = Region( "1", 1, 2 );
    auto region2 = Region( "1", 2, 4 );
    auto region3 = Region( "1", 3, 3 );
    auto region4 = Region( "2", 2, 3 );
    auto refSequenceChr1 =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'A' ) );
    auto refSequenceChr2 =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "2", 0, 5 ), std::string( 5, 'A' ) );

    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequenceChr1, region1, "C" );
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequenceChr1, region2, "TC" );
    echidna::variant::varPtr_t variant3 = std::make_shared< Variant >( refSequenceChr1, region3, "CTA" );
    echidna::variant::varPtr_t variant4 = std::make_shared< Variant >( refSequenceChr2, region4, "C" );

    Call call1( variant1, variant1->interval(), 10, 1, {{Call::REF, Call::VAR}} );
    Call call2( variant2, variant2->interval(), 10, 1, {{Call::REF, Call::VAR}} );
    Call call3( variant3, variant3->interval(), 10, 1, {{Call::REF, Call::VAR}} );
    Call call4( variant4, variant4->interval(), 10, 1, {{Call::REF, Call::VAR}} );

    const auto callMerger = echidna::caller::lvcMerger();

    auto result_no_overlap = callMerger.getOverlappingLvcCallIndices( call1, {call2, call3, call4} );
    BOOST_CHECK_EQUAL( result_no_overlap.size(), 0 );

    auto result_overlap = callMerger.getOverlappingLvcCallIndices( call2, {call1, call3, call4} );
    BOOST_CHECK_EQUAL( result_overlap.size(), 1 );
    BOOST_CHECK_EQUAL( result_overlap[0], 1 );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( shouldMergeHomozogousLargeVariantWithOtherVariants )
{
    auto region1 = Region( "1", 1, 2 );
    auto region2 = Region( "1", 2, 80 );
    auto region3 = Region( "1", 30, 31 );
    auto refSequenceChr1 =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 100 ), std::string( 100, 'A' ) );

    echidna::variant::varPtr_t lvcSNP = std::make_shared< Variant >( refSequenceChr1, region1, "C" );
    echidna::variant::varPtr_t lvcLargeDeletion = std::make_shared< Variant >( refSequenceChr1, region2, "" );
    echidna::variant::varPtr_t variant = std::make_shared< Variant >( refSequenceChr1, region3, "C" );

    Call call1( lvcSNP, lvcSNP->interval(), 10, 1, {{Call::REF, Call::VAR}} );
    Call call2( lvcLargeDeletion, lvcLargeDeletion->interval(), 10, 1, {{Call::VAR, Call::VAR}} );
    Call call3( variant, variant->interval(), 10, 1, {{}} );
    BOOST_CHECK_EQUAL( call3.samples[0].genotypeCalls.size(), 0 );

    callVector_t lvcCalls = {call1, call2};
    callVector_t overlappingCalls = {call3};

    const auto callMerger = echidna::caller::lvcMerger();

    callMerger.mergeAndCorrectGenotypes( overlappingCalls, lvcCalls, {0}, {2} );

    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls[0], Call::UNKNOWN );
    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls[1], Call::UNKNOWN );
}

BOOST_AUTO_TEST_CASE( shouldMergeHetLargeVariantOnSecondStrandWithOtherVariants )
{
    auto region1 = Region( "1", 1, 2 );
    auto region2 = Region( "1", 2, 80 );
    auto region3 = Region( "1", 30, 31 );
    auto refSequenceChr1 =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 100 ), std::string( 100, 'A' ) );

    echidna::variant::varPtr_t lvcSNP = std::make_shared< Variant >( refSequenceChr1, region1, "C" );
    echidna::variant::varPtr_t lvcLargeDeletion = std::make_shared< Variant >( refSequenceChr1, region2, "" );
    echidna::variant::varPtr_t variant = std::make_shared< Variant >( refSequenceChr1, region3, "C" );

    Call call1( lvcSNP, lvcSNP->interval(), 10, 1, {{Call::REF, Call::VAR}} );
    Call call2( lvcLargeDeletion, lvcLargeDeletion->interval(), 10, 1, {{Call::REF, Call::VAR}} );
    Call call3( variant, variant->interval(), 10, 1, {{Call::VAR}} );  // variant for ploidy=1

    callVector_t lvcCalls = {call1, call2};
    callVector_t overlappingCalls = {call3};

    const auto callMerger = echidna::caller::lvcMerger();

    callMerger.mergeAndCorrectGenotypes( overlappingCalls, lvcCalls, {1}, {2} );

    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls[0], Call::VAR );
    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls[1], Call::UNKNOWN );

    BOOST_CHECK_EQUAL( lvcCalls[1].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( lvcCalls[1].samples[0].genotypeCalls[0], Call::UNKNOWN );
    BOOST_CHECK_EQUAL( lvcCalls[1].samples[0].genotypeCalls[1], Call::VAR );
}

BOOST_AUTO_TEST_CASE( shouldMergeHetLargeVariantOnFirstStrandWithOtherVariants )
{
    auto region1 = Region( "1", 1, 2 );
    auto region2 = Region( "1", 2, 80 );
    auto region3 = Region( "1", 30, 31 );
    auto refSequenceChr1 =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 100 ), std::string( 100, 'A' ) );

    echidna::variant::varPtr_t lvcSNP = std::make_shared< Variant >( refSequenceChr1, region1, "C" );
    echidna::variant::varPtr_t lvcLargeDeletion = std::make_shared< Variant >( refSequenceChr1, region2, "" );
    echidna::variant::varPtr_t variant = std::make_shared< Variant >( refSequenceChr1, region3, "C" );

    Call call1( lvcSNP, lvcSNP->interval(), 10, 1, {{Call::REF, Call::VAR}} );
    Call call2( lvcLargeDeletion, lvcLargeDeletion->interval(), 10, 1, {{Call::VAR, Call::REF}} );
    Call call3( variant, variant->interval(), 10, 1, {{Call::VAR}} );  // variant for ploidy=1

    callVector_t lvcCalls = {call1, call2};
    callVector_t overlappingCalls = {call3};

    const auto callMerger = echidna::caller::lvcMerger();

    callMerger.mergeAndCorrectGenotypes( overlappingCalls, lvcCalls, {1}, {2} );

    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls[0], Call::UNKNOWN );
    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls[1], Call::VAR );

    BOOST_CHECK_EQUAL( lvcCalls[1].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( lvcCalls[1].samples[0].genotypeCalls[0], Call::VAR );
    BOOST_CHECK_EQUAL( lvcCalls[1].samples[0].genotypeCalls[1], Call::UNKNOWN );
}

BOOST_AUTO_TEST_CASE( shouldMergeHomozogousLargeVariantWithTwoOverlappingLargeVariants )
{
    auto region1 = Region( "1", 1, 2 );
    auto lvRegion1 = Region( "1", 2, 80 );
    auto lvRegion2 = Region( "1", 10, 90 );
    auto region2 = Region( "1", 30, 31 );
    auto refSequenceChr1 =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 100 ), std::string( 100, 'A' ) );

    echidna::variant::varPtr_t lvcSNP = std::make_shared< Variant >( refSequenceChr1, region1, "C" );
    echidna::variant::varPtr_t lvcLargeDeletion1 = std::make_shared< Variant >( refSequenceChr1, lvRegion1, "" );
    echidna::variant::varPtr_t lvcLargeDeletion2 = std::make_shared< Variant >( refSequenceChr1, lvRegion2, "" );
    echidna::variant::varPtr_t variant = std::make_shared< Variant >( refSequenceChr1, region2, "C" );

    Call call1( lvcSNP, lvcSNP->interval(), 10, 1, {{Call::REF, Call::VAR}} );
    Call call2( lvcLargeDeletion1, lvcLargeDeletion1->interval(), 10, 1, {{Call::UNKNOWN, Call::UNKNOWN}} );
    Call call3( lvcLargeDeletion2, lvcLargeDeletion2->interval(), 10, 1, {{Call::VAR, Call::REF}} );
    Call call4( variant, variant->interval(), 10, 1, {{Call::VAR}} );
    BOOST_CHECK_EQUAL( call4.samples[0].genotypeCalls.size(), 1 );

    callVector_t lvcCalls = {call1, call2, call3};
    callVector_t overlappingCalls = {call4};

    const auto callMerger = echidna::caller::lvcMerger();

    callMerger.mergeAndCorrectGenotypes( overlappingCalls, lvcCalls, {1}, {2} );

    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls[0], Call::UNKNOWN );
    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls[1], Call::VAR );

    // first large variant: not called
    BOOST_CHECK_EQUAL( lvcCalls[1].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( lvcCalls[1].samples[0].genotypeCalls[0], Call::UNKNOWN );
    BOOST_CHECK_EQUAL( lvcCalls[1].samples[0].genotypeCalls[1], Call::UNKNOWN );

    // second large variant: het
    BOOST_CHECK_EQUAL( lvcCalls[2].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( lvcCalls[2].samples[0].genotypeCalls[0], Call::VAR );
    BOOST_CHECK_EQUAL( lvcCalls[2].samples[0].genotypeCalls[1], Call::UNKNOWN );
}

BOOST_AUTO_TEST_CASE( shouldThrowForInconsistentGenotypes )
{
    auto region2 = Region( "1", 2, 80 );
    auto region3 = Region( "1", 30, 31 );
    auto refSequenceChr1 =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 100 ), std::string( 100, 'A' ) );

    echidna::variant::varPtr_t lvcLargeDeletion = std::make_shared< Variant >( refSequenceChr1, region2, "" );
    echidna::variant::varPtr_t variant = std::make_shared< Variant >( refSequenceChr1, region3, "C" );

    Call call2( lvcLargeDeletion, lvcLargeDeletion->interval(), 10, 1, {{Call::VAR, Call::VAR}} );
    Call call3( variant, variant->interval(), 10, 1, {{Call::VAR}} );  // variant for ploidy=1

    callVector_t lvcCalls = {call2};
    callVector_t overlappingCalls = {call3};

    const auto callMerger = echidna::caller::lvcMerger();

    BOOST_CHECK_THROW( callMerger.mergeAndCorrectGenotypes( overlappingCalls, lvcCalls, {1}, {2} ),
                       echidna::utils::echidna_exception );
}

BOOST_AUTO_TEST_CASE( shouldMergeHetLargeVariantForTwoSamplesWhereSecondSampleIsHom )
{
    // This is testing the merging of a large variant cluster for two samples and 3 variants.
    // sample 1 - strand1: large deletion
    // sample 1 - strand2: snp, variant overlapping large deletion
    // sample 2 - strand1: large deletion
    // sample 2 - strand2: large deletion

    auto nSamples = 2;
    auto region1 = Region( "1", 1, 2 );
    auto region2 = Region( "1", 2, 80 );
    auto region3 = Region( "1", 30, 31 );
    auto refSequenceChr1 =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 100 ), std::string( 100, 'A' ) );

    echidna::variant::varPtr_t lvcSNP = std::make_shared< Variant >( refSequenceChr1, region1, "C" );
    echidna::variant::varPtr_t lvcLargeDeletion = std::make_shared< Variant >( refSequenceChr1, region2, "" );
    echidna::variant::varPtr_t variant = std::make_shared< Variant >( refSequenceChr1, region3, "C" );

    Call call1( lvcSNP, lvcSNP->interval(), 10, nSamples, {{Call::REF, Call::VAR}, {Call::REF, Call::REF}} );
    Call call2( lvcLargeDeletion, lvcLargeDeletion->interval(), 10, nSamples,
                {{Call::VAR, Call::REF}, {Call::VAR, Call::VAR}} );
    Call call3( variant, variant->interval(), 10, nSamples, {{Call::VAR}, {}} );  // variant with reduced ploidy

    callVector_t lvcCalls = {call1, call2};
    callVector_t overlappingCalls = {call3};

    const auto callMerger = echidna::caller::lvcMerger();

    callMerger.mergeAndCorrectGenotypes( overlappingCalls, lvcCalls, {1, 0}, {2, 2} );

    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls[0], Call::UNKNOWN );
    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls[1], Call::VAR );

    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[1].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[1].genotypeCalls[0], Call::UNKNOWN );
    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[1].genotypeCalls[1], Call::UNKNOWN );

    BOOST_CHECK_EQUAL( lvcCalls[1].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( lvcCalls[1].samples[0].genotypeCalls[0], Call::VAR );
    BOOST_CHECK_EQUAL( lvcCalls[1].samples[0].genotypeCalls[1], Call::UNKNOWN );

    BOOST_CHECK_EQUAL( lvcCalls[1].samples[1].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( lvcCalls[1].samples[1].genotypeCalls[0], Call::VAR );
    BOOST_CHECK_EQUAL( lvcCalls[1].samples[1].genotypeCalls[1], Call::VAR );
}

BOOST_AUTO_TEST_CASE( shouldMergeHetLargeVariantForTwoSamplesWhereFirstSampleIsHom )
{
    // This is testing the merging of a large variant cluster for two samples and 3 variants.
    // sample 1 - strand1: large deletion
    // sample 1 - strand2: snp, large deletion
    // sample 2 - strand1: snp, variant overlapping large deletion
    // sample 2 - strand2: large deletion

    auto nSamples = 2;
    auto region1 = Region( "1", 1, 2 );
    auto region2 = Region( "1", 2, 80 );
    auto region3 = Region( "1", 30, 31 );
    auto refSequenceChr1 =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 100 ), std::string( 100, 'A' ) );

    echidna::variant::varPtr_t lvcSNP = std::make_shared< Variant >( refSequenceChr1, region1, "C" );
    echidna::variant::varPtr_t lvcLargeDeletion = std::make_shared< Variant >( refSequenceChr1, region2, "" );
    echidna::variant::varPtr_t variant = std::make_shared< Variant >( refSequenceChr1, region3, "C" );

    Call call1( lvcSNP, lvcSNP->interval(), 10, nSamples, {{Call::REF, Call::VAR}, {Call::VAR, Call::REF}} );
    Call call2( lvcLargeDeletion, lvcLargeDeletion->interval(), 10, nSamples,
                {{Call::VAR, Call::VAR}, {Call::REF, Call::VAR}} );
    Call call3( variant, variant->interval(), 10, nSamples, {{}, {Call::VAR}} );

    callVector_t lvcCalls = {call1, call2};
    callVector_t overlappingCalls = {call3};

    const auto callMerger = echidna::caller::lvcMerger();

    callMerger.mergeAndCorrectGenotypes( overlappingCalls, lvcCalls, {0, 1}, {2, 2} );

    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls[0], Call::UNKNOWN );
    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[0].genotypeCalls[1], Call::UNKNOWN );

    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[1].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[1].genotypeCalls[0], Call::VAR );
    BOOST_CHECK_EQUAL( overlappingCalls[0].samples[1].genotypeCalls[1], Call::UNKNOWN );

    BOOST_CHECK_EQUAL( lvcCalls[1].samples[0].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( lvcCalls[1].samples[0].genotypeCalls[0], Call::VAR );
    BOOST_CHECK_EQUAL( lvcCalls[1].samples[0].genotypeCalls[1], Call::VAR );

    BOOST_CHECK_EQUAL( lvcCalls[1].samples[1].genotypeCalls.size(), 2 );
    BOOST_CHECK_EQUAL( lvcCalls[1].samples[1].genotypeCalls[0], Call::UNKNOWN );
    BOOST_CHECK_EQUAL( lvcCalls[1].samples[1].genotypeCalls[1], Call::VAR );
}