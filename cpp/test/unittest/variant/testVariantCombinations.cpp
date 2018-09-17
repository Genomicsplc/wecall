#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "variant/type/variant.hpp"
#include "variant/variantCombinations.hpp"
#include "variant/clustering.hpp"
#include "io/read.hpp"
#include "alignment/cigar.hpp"
#include "caller/region.hpp"

using namespace echidna::variant;
using namespace echidna::alignment;
using namespace echidna::caller;
using namespace echidna::io;
using echidna::utils::ReferenceSequence;

readPtr_t makeRead( int64_t startPos, int64_t endPos )
{
    const auto length = int64_to_sizet( endPos - startPos );
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", startPos - 10, endPos + 10 ),
                                                                          std::string( length + 20, 'A' ) );

    return std::make_shared< Read >( echidna::utils::BasePairSequence( length, 'T' ), std::string( length, 'Q' ), "",
                                     Cigar( std::to_string( length ) + std::string( "M" ) ), 0, startPos, 0, 0, 0, 0, 0,
                                     referenceSequence );
}

// --- compute variant combinations ------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testShouldGetVariantComboForReference )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 15 ), "AAAAA" );
    VariantCombinations combos( 1, 40 );

    BOOST_CHECK_EQUAL( combos.nVariantCombinations(), 0 );
    combos.computeVariantCombinations( {}, 10, referenceSequence );
    BOOST_REQUIRE_EQUAL( combos.nVariantCombinations(), 1 );
    BOOST_CHECK_EQUAL( combos.variantCombinations().size(), 1 );
    BOOST_CHECK_EQUAL( combos.variantCombinations()[0].size(), 0 );
}

BOOST_AUTO_TEST_CASE( testShouldReduceNumberOfCombinationsIfVariantIsCorrelatedToPreviousVariant )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 15 ), "AAAAA" );
    VariantCluster cluster;

    auto read1 = std::make_shared< Read >( "ATG", "QQQ", "", Cigar( "3M" ), 0, 11, 0, 0, 0, 0, 0, referenceSequence );
    auto read2 = std::make_shared< Read >( "ATG", "QQQ", "", Cigar( "3M" ), 0, 11, 0, 0, 0, 0, 0, referenceSequence );

    const auto snp1 = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "T" );
    cluster.push_back( snp1, Region( "1", 11, 12 ) );

    const auto snp2 = std::make_shared< Variant >( referenceSequence, Region( "1", 12, 13 ), "T" );
    snp2->addRead( read1 );
    snp2->addRead( read2 );
    cluster.push_back( snp2, Region( "1", 12, 13 ) );

    const auto snp3 = std::make_shared< Variant >( referenceSequence, Region( "1", 13, 14 ), "T" );
    snp3->addRead( read1 );
    snp3->addRead( read2 );
    cluster.push_back( snp3, Region( "1", 13, 14 ) );

    const auto maxCombinations = 100;

    VariantCombinations combos( 1, 40 );
    BOOST_REQUIRE_EQUAL( combos.nVariantCombinations(), 0 );

    combos.computeVariantCombinations( cluster.variants(), maxCombinations, referenceSequence );

    BOOST_REQUIRE_EQUAL( combos.nVariantCombinations(), 4 );
    BOOST_REQUIRE_EQUAL( combos.variantCombinations().size(), 4 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[0].size(), 1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[0][0], *snp1 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[1].size(), 3 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[1][0], *snp1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[1][1], *snp2 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[1][2], *snp3 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[2].size(), 2 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[2][0], *snp2 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[2][1], *snp3 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[3].size(), 0 );
}

BOOST_AUTO_TEST_CASE( testShouldOnlyGetIsolatedSNPCombosIfVariantIsUncorrelatedToPreviousVariant )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 15 ), "AAAAA" );
    VariantCluster cluster;

    auto read1 = std::make_shared< Read >( "ATG", "QQQ", "", Cigar( "3M" ), 0, 11, 0, 0, 0, 0, 0, referenceSequence );
    auto read2 = std::make_shared< Read >( "ATG", "QQQ", "", Cigar( "3M" ), 0, 11, 0, 0, 0, 0, 0, referenceSequence );

    const auto snp1 = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "T" );
    cluster.push_back( snp1, Region( "1", 11, 12 ) );

    const auto snp2 = std::make_shared< Variant >( referenceSequence, Region( "1", 12, 13 ), "T" );
    snp2->addRead( read1 );
    cluster.push_back( snp2, Region( "1", 12, 13 ) );

    const auto snp3 = std::make_shared< Variant >( referenceSequence, Region( "1", 13, 14 ), "T" );
    snp3->addRead( read2 );
    cluster.push_back( snp3, Region( "1", 13, 14 ) );

    const auto maxCombinations = 100;

    VariantCombinations combos( 1, 40 );
    BOOST_CHECK_EQUAL( combos.nVariantCombinations(), 0 );

    combos.computeVariantCombinations( cluster.variants(), maxCombinations, referenceSequence );
    BOOST_REQUIRE_EQUAL( combos.nVariantCombinations(), 6 );
    BOOST_CHECK_EQUAL( combos.variantCombinations().size(), 6 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[0].size(), 1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[0][0], *snp1 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[1].size(), 2 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[1][0], *snp1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[1][1], *snp2 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[2].size(), 1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[2][0], *snp2 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[3].size(), 2 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[3][0], *snp1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[3][1], *snp3 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[4].size(), 1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[4][0], *snp3 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[5].size(), 0 );
}

BOOST_AUTO_TEST_CASE( testShouldOnlyGetCorrectCombinationsIfFirstVariantImpliesSecondVariant )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 15 ), "AAAAA" );
    VariantCluster cluster;

    auto read1 = std::make_shared< Read >( "ATG", "QQQ", "", Cigar( "3M" ), 0, 11, 0, 0, 0, 0, 0, referenceSequence );
    auto read2 = std::make_shared< Read >( "ATG", "QQQ", "", Cigar( "3M" ), 0, 11, 0, 0, 0, 0, 0, referenceSequence );
    std::cout << "\n";

    const auto snp1 = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "T" );
    cluster.push_back( snp1, Region( "1", 11, 12 ) );

    const auto snp2 = std::make_shared< Variant >( referenceSequence, Region( "1", 12, 13 ), "T" );
    snp2->addRead( read1 );
    cluster.push_back( snp2, Region( "1", 12, 13 ) );

    const auto snp3 = std::make_shared< Variant >( referenceSequence, Region( "1", 13, 14 ), "T" );
    snp3->addRead( read1 );
    snp3->addRead( read2 );
    cluster.push_back( snp3, Region( "1", 13, 14 ) );

    const auto maxCombinations = 100;

    VariantCombinations combos( 1, 40 );
    BOOST_CHECK_EQUAL( combos.nVariantCombinations(), 0 );

    combos.computeVariantCombinations( cluster.variants(), maxCombinations, referenceSequence );

    BOOST_REQUIRE_EQUAL( combos.nVariantCombinations(), 6 );
    BOOST_CHECK_EQUAL( combos.variantCombinations().size(), 6 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[0].size(), 1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[0][0], *snp1 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[1].size(), 3 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[1][0], *snp1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[1][1], *snp2 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[1][2], *snp3 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[2].size(), 2 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[2][0], *snp2 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[2][1], *snp3 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[3].size(), 2 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[3][0], *snp1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[3][1], *snp3 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[4].size(), 1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[4][0], *snp3 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[5].size(), 0 );
}

BOOST_AUTO_TEST_CASE( testShouldOnlyGetCorrectCombinationsIfSecondVariantImpliesFirstVariant )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 15 ), "AAAAA" );
    VariantCluster cluster;

    auto read1 = std::make_shared< Read >( "ATG", "QQQ", "", Cigar( "3M" ), 0, 11, 0, 0, 0, 0, 0, referenceSequence );
    auto read2 = std::make_shared< Read >( "ATG", "QQQ", "", Cigar( "3M" ), 0, 11, 0, 0, 0, 0, 0, referenceSequence );
    std::cout << "\n";

    const auto snp1 = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "T" );
    cluster.push_back( snp1, Region( "1", 11, 12 ) );

    const auto snp2 = std::make_shared< Variant >( referenceSequence, Region( "1", 12, 13 ), "T" );
    snp2->addRead( read1 );
    snp2->addRead( read2 );
    cluster.push_back( snp2, Region( "1", 12, 13 ) );

    const auto snp3 = std::make_shared< Variant >( referenceSequence, Region( "1", 13, 14 ), "T" );
    snp3->addRead( read1 );
    cluster.push_back( snp3, Region( "1", 13, 14 ) );

    const auto maxCombinations = 100;

    VariantCombinations combos( 1, 40 );
    BOOST_CHECK_EQUAL( combos.nVariantCombinations(), 0 );

    combos.computeVariantCombinations( cluster.variants(), maxCombinations, referenceSequence );

    BOOST_REQUIRE_EQUAL( combos.nVariantCombinations(), 7 );
    BOOST_CHECK_EQUAL( combos.variantCombinations().size(), 7 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[0].size(), 1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[0][0], *snp1 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[1].size(), 2 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[1][0], *snp1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[1][1], *snp2 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[2].size(), 1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[2][0], *snp2 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[3].size(), 2 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[3][0], *snp1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[3][1], *snp3 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[4].size(), 3 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[4][0], *snp1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[4][1], *snp2 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[4][2], *snp3 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[5].size(), 2 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[5][0], *snp2 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[5][1], *snp3 );

    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[6].size(), 0 );
}

BOOST_AUTO_TEST_CASE( testShouldLimitNumberOfCombinationsWhenPassedLimit )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 15 ), "AAAAA" );
    VariantCluster cluster;

    const auto snp1 = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "T" );
    cluster.push_back( snp1, Region( "1", 11, 12 ) );

    const auto snp2 = std::make_shared< Variant >( referenceSequence, Region( "1", 12, 13 ), "T" );
    cluster.push_back( snp2, Region( "1", 12, 13 ) );

    const auto maxCombinations = 2;

    VariantCombinations combos( 1, 40 );

    BOOST_CHECK_EQUAL( combos.nVariantCombinations(), 0 );
    BOOST_CHECK( not combos.allCombinationsComputed() );

    combos.computeVariantCombinations( cluster.variants(), maxCombinations, referenceSequence );

    BOOST_CHECK( not combos.allCombinationsComputed() );
    BOOST_CHECK_EQUAL( combos.nVariantCombinations(), 0 );
    BOOST_CHECK_EQUAL( combos.variantCombinations().size(), 0 );

    combos.computeVariantCombinations( cluster.variants(), 100, referenceSequence );
    BOOST_CHECK( combos.allCombinationsComputed() );
    BOOST_REQUIRE_EQUAL( combos.nVariantCombinations(), 4 );
    BOOST_REQUIRE_EQUAL( combos.variantCombinations().size(), 4 );
    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[0].size(), 1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[0][0], *snp1 );
    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[1].size(), 2 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[1][0], *snp1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[1][1], *snp2 );
    BOOST_REQUIRE_EQUAL( combos.variantCombinations()[2].size(), 1 );
    BOOST_CHECK_EQUAL( *combos.variantCombinations()[2][0], *snp2 );
    BOOST_CHECK_EQUAL( combos.variantCombinations()[3].size(), 0 );
}

// --- filter variant combinations -----------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testShouldFilterInvalidCombinationsForAlwaysTogetherCase )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 15 ), "AAAAA" );
    const auto snp1 = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "T" );
    const auto snp2 = std::make_shared< Variant >( referenceSequence, Region( "1", 12, 13 ), "T" );
    const auto snp3 = std::make_shared< Variant >( referenceSequence, Region( "1", 13, 14 ), "T" );

    std::vector< std::vector< varPtr_t > > variantCombinations;
    variantCombinations.push_back( {snp3} );              // invalid
    variantCombinations.push_back( {snp1, snp3} );        // invalid
    variantCombinations.push_back( {snp1, snp2, snp3} );  // valid

    std::vector< varPtr_t > alwaysTogetherVariants;
    alwaysTogetherVariants.push_back( snp2 );

    std::vector< varPtr_t > neverTogetherVariants;

    VariantCombinations combos( 1, 40 );
    auto validVariantCombinations =
        combos.filterVariantCombinations( variantCombinations, alwaysTogetherVariants, neverTogetherVariants, snp3 );
    BOOST_REQUIRE_EQUAL( validVariantCombinations.size(), 1 );
    BOOST_REQUIRE_EQUAL( validVariantCombinations[0].size(), 3 );
    BOOST_CHECK_EQUAL( *validVariantCombinations[0][0], *snp1 );
    BOOST_CHECK_EQUAL( *validVariantCombinations[0][1], *snp2 );
    BOOST_CHECK_EQUAL( *validVariantCombinations[0][2], *snp3 );
}

BOOST_AUTO_TEST_CASE( testShouldFilterInvalidCombinationsForNeverTogetherCase )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 15 ), "AAAAA" );
    const auto snp1 = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "T" );
    const auto snp2 = std::make_shared< Variant >( referenceSequence, Region( "1", 12, 13 ), "T" );
    const auto snp3 = std::make_shared< Variant >( referenceSequence, Region( "1", 13, 14 ), "T" );

    std::vector< std::vector< varPtr_t > > variantCombinations;
    variantCombinations.push_back( {snp1} );              // valid
    variantCombinations.push_back( {snp3} );              // valid
    variantCombinations.push_back( {snp2, snp3} );        // invalid
    variantCombinations.push_back( {snp1, snp2, snp3} );  // invalid

    std::vector< varPtr_t > alwaysTogetherVariants;

    std::vector< varPtr_t > neverTogetherVariants;
    neverTogetherVariants.push_back( snp2 );

    VariantCombinations combos( 1, 40 );
    auto validVariantCombinations =
        combos.filterVariantCombinations( variantCombinations, alwaysTogetherVariants, neverTogetherVariants, snp3 );

    BOOST_REQUIRE_EQUAL( validVariantCombinations.size(), 2 );

    BOOST_REQUIRE_EQUAL( validVariantCombinations[0].size(), 1 );
    BOOST_CHECK_EQUAL( *validVariantCombinations[0][0], *snp1 );

    BOOST_REQUIRE_EQUAL( validVariantCombinations[1].size(), 1 );
    BOOST_CHECK_EQUAL( *validVariantCombinations[1][0], *snp3 );
}

BOOST_AUTO_TEST_CASE( testShouldFilterInvalidCombinationsForMixedCase )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 15 ), "AAAAA" );
    const auto snp1 = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "T" );
    const auto snp2 = std::make_shared< Variant >( referenceSequence, Region( "1", 12, 13 ), "T" );
    const auto snp3 = std::make_shared< Variant >( referenceSequence, Region( "1", 13, 14 ), "T" );
    const auto snp4 = std::make_shared< Variant >( referenceSequence, Region( "1", 14, 15 ), "T" );

    std::vector< std::vector< varPtr_t > > variantCombinations;
    variantCombinations.push_back( {snp1, snp4} );        // invalid
    variantCombinations.push_back( {snp1, snp3, snp4} );  // invalid
    variantCombinations.push_back( {snp4} );              // invalid
    variantCombinations.push_back( {snp3} );              // valid
    variantCombinations.push_back( {snp2, snp4} );        // valid
    variantCombinations.push_back( {snp2, snp3, snp4} );  // valid

    std::vector< varPtr_t > alwaysTogetherVariants;
    alwaysTogetherVariants.push_back( snp2 );

    std::vector< varPtr_t > neverTogetherVariants;
    neverTogetherVariants.push_back( snp1 );

    VariantCombinations combos( 1, 40 );
    auto validVariantCombinations =
        combos.filterVariantCombinations( variantCombinations, alwaysTogetherVariants, neverTogetherVariants, snp4 );
    BOOST_REQUIRE_EQUAL( validVariantCombinations.size(), 3 );

    BOOST_REQUIRE_EQUAL( validVariantCombinations[0].size(), 1 );
    BOOST_CHECK_EQUAL( *validVariantCombinations[0][0], *snp3 );

    BOOST_REQUIRE_EQUAL( validVariantCombinations[1].size(), 2 );
    BOOST_CHECK_EQUAL( *validVariantCombinations[1][0], *snp2 );
    BOOST_CHECK_EQUAL( *validVariantCombinations[1][1], *snp4 );

    BOOST_REQUIRE_EQUAL( validVariantCombinations[2].size(), 3 );
    BOOST_CHECK_EQUAL( *validVariantCombinations[2][0], *snp2 );
    BOOST_CHECK_EQUAL( *validVariantCombinations[2][1], *snp3 );
    BOOST_CHECK_EQUAL( *validVariantCombinations[2][2], *snp4 );
}

// --- get state for variants ----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( test_should_return_uncertain_for_variants_with_no_reads )
{
    auto reference = std::make_shared< ReferenceSequence >( Region( "1", 0, 2 ), "AA" );
    auto first = std::make_shared< Variant >( reference, Region( "1", 0, 1 ), "T" );
    auto second = std::make_shared< Variant >( reference, Region( "1", 1, 2 ), "T" );

    BOOST_CHECK_EQUAL( VariantPairCombinationState::UNCERTAIN, VariantCombinations( 1, 40 ).getState( first, second ) );
}

BOOST_AUTO_TEST_CASE( test_should_return_state_always_together_if_all_reads_support_both_variants )
{
    auto reference = std::make_shared< ReferenceSequence >( Region( "1", 0, 2 ), "AA" );
    auto first = std::make_shared< Variant >( reference, Region( "1", 0, 1 ), "T" );
    auto second = std::make_shared< Variant >( reference, Region( "1", 1, 2 ), "T" );

    auto read1 = makeRead( 0, 3 );
    auto read2 = makeRead( 0, 3 );

    first->addRead( read2 );
    first->addRead( read1 );

    second->addRead( read1 );
    second->addRead( read2 );

    BOOST_CHECK_EQUAL( VariantPairCombinationState::ALWAYS_TOGETHER,
                       VariantCombinations( 1, 40 ).getState( first, second ) );
}

BOOST_AUTO_TEST_CASE( test_should_return_state_uncertain_if_not_enough_reads_support_the_claim )
{
    auto reference = std::make_shared< ReferenceSequence >( Region( "1", 0, 2 ), "AA" );
    auto first = std::make_shared< Variant >( reference, Region( "1", 0, 1 ), "T" );
    auto second = std::make_shared< Variant >( reference, Region( "1", 1, 2 ), "T" );

    auto read1 = makeRead( 0, 3 );
    auto read2 = makeRead( 0, 3 );

    first->addRead( read2 );
    first->addRead( read1 );

    second->addRead( read1 );
    second->addRead( read2 );

    std::size_t nNeededReads = 3;
    BOOST_CHECK_EQUAL( VariantPairCombinationState::UNCERTAIN,
                       VariantCombinations( nNeededReads, 40 ).getState( first, second ) );
}

BOOST_AUTO_TEST_CASE( test_should_return_state_never_together_if_all_reads_support_different_variants )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 2 ), "AA" );
    auto first = std::make_shared< Variant >( reference, Region( "1", startPos, startPos + 1 ), "T" );
    auto second = std::make_shared< Variant >( reference, Region( "1", startPos + 1, startPos + 2 ), "T" );

    auto read1 = makeRead( 0, 3 );
    auto read2 = makeRead( 0, 3 );

    // Maybe add more data!
    first->addRead( read1 );
    second->addRead( read2 );

    BOOST_CHECK_EQUAL( VariantPairCombinationState::NEVER_TOGETHER,
                       VariantCombinations( 1, 40 ).getState( first, second ) );
}

BOOST_AUTO_TEST_CASE( test_should_return_state_first_implies_second_if_all_reads_support_first_implies_second_variant )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 2 ), "AA" );
    auto first = std::make_shared< Variant >( reference, Region( "1", startPos, startPos + 1 ), "T" );
    auto second = std::make_shared< Variant >( reference, Region( "1", startPos + 1, startPos + 2 ), "T" );

    auto read1 = makeRead( 0, 3 );
    auto read2 = makeRead( 0, 3 );

    // Maybe add more data!
    first->addRead( read1 );
    first->addRead( read2 );

    second->addRead( read2 );

    BOOST_CHECK_EQUAL( VariantPairCombinationState::SECOND_IMPLIES_FIRST,
                       VariantCombinations( 1, 40 ).getState( first, second ) );
}

BOOST_AUTO_TEST_CASE( test_should_return_state_second_implies_first_if_all_reads_support_second_implies_first_variant )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 2 ), "AA" );
    auto first = std::make_shared< Variant >( reference, Region( "1", startPos, startPos + 1 ), "T" );
    auto second = std::make_shared< Variant >( reference, Region( "1", startPos + 1, startPos + 2 ), "T" );

    auto read1 = makeRead( 0, 3 );
    auto read2 = makeRead( 0, 3 );

    // Maybe add more data!
    first->addRead( read1 );

    second->addRead( read1 );
    second->addRead( read2 );

    BOOST_CHECK_EQUAL( VariantPairCombinationState::FIRST_IMPLIES_SECOND,
                       VariantCombinations( 1, 40 ).getState( first, second ) );
}

BOOST_AUTO_TEST_CASE( test_should_return_uncertain_for_mixed_read_situation )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 2 ), "AA" );
    auto first = std::make_shared< Variant >( reference, Region( "1", startPos, startPos + 1 ), "T" );
    auto second = std::make_shared< Variant >( reference, Region( "1", startPos + 1, startPos + 2 ), "T" );

    auto read1 = makeRead( startPos, 3 );
    auto read2 = makeRead( startPos, 3 );
    auto read3 = makeRead( startPos, 3 );

    // Maybe add more data!
    first->addRead( read1 );
    first->addRead( read3 );

    second->addRead( read1 );
    second->addRead( read2 );

    BOOST_CHECK_EQUAL( VariantPairCombinationState::UNCERTAIN, VariantCombinations( 1, 40 ).getState( first, second ) );
}

BOOST_AUTO_TEST_CASE( test_should_return_uncertain_for_one_variant_without_overlapping_reads )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 2 ), "AA" );
    auto first = std::make_shared< Variant >( reference, Region( "1", startPos, startPos + 1 ), "T" );
    auto second = std::make_shared< Variant >( reference, Region( "1", startPos + 1, startPos + 2 ), "T" );

    auto read1 = makeRead( startPos, 3 );

    second->addRead( read1 );

    BOOST_CHECK_EQUAL( VariantPairCombinationState::UNCERTAIN, VariantCombinations( 1, 40 ).getState( first, second ) );
}

BOOST_AUTO_TEST_CASE( test_should_return_always_together_for_mixed_read_situation_considering_read_regions )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 2 ), "AA" );

    auto first = std::make_shared< Variant >( reference, Region( "1", startPos, startPos + 1 ), "T" );
    auto second = std::make_shared< Variant >( reference, Region( "1", startPos + 1, startPos + 2 ), "T" );

    auto read1 = makeRead( startPos, startPos + 1 );
    auto read2 = makeRead( startPos + 1, startPos + 2 );
    auto read3 = makeRead( startPos, startPos + 3 );

    // Maybe add more data!
    first->addRead( read1 );
    first->addRead( read3 );

    second->addRead( read2 );
    second->addRead( read3 );

    VariantPairCombinationState state = VariantCombinations( 1, 40 ).getState( first, second );
    BOOST_CHECK_EQUAL( VariantPairCombinationState::ALWAYS_TOGETHER, state );
}

// --- filter reads -------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( test_should_apply_read_filter_to_empty_case )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 2 ), "AA" );

    auto first = std::make_shared< Variant >( reference, Region( "1", startPos, startPos + 1 ), "T" );
    auto second = std::make_shared< Variant >( reference, Region( "1", startPos + 1, startPos + 2 ), "T" );

    std::vector< readPtr_t > firstReads;
    std::vector< readPtr_t > secondReads;

    VariantCombinations( 1, 40 ).getOverlappingReads( first, second, firstReads, secondReads );
    BOOST_CHECK( firstReads.empty() );
    BOOST_CHECK( secondReads.empty() );
}

BOOST_AUTO_TEST_CASE( test_should_return_no_reads_if_first_variant_has_no_reads )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 2 ), "AA" );
    auto read1 = makeRead( startPos, startPos + 1 );

    auto first = std::make_shared< Variant >( reference, Region( "1", startPos, startPos + 1 ), "T" );
    first->addRead( read1 );

    auto second = std::make_shared< Variant >( reference, Region( "1", startPos + 1, startPos + 2 ), "T" );

    std::vector< readPtr_t > firstReads;
    std::vector< readPtr_t > secondReads;
    VariantCombinations( 1, 40 ).getOverlappingReads( first, second, firstReads, secondReads );
    BOOST_CHECK( firstReads.empty() );
    BOOST_CHECK( secondReads.empty() );
}

BOOST_AUTO_TEST_CASE( test_should_return_no_reads_if_second_variant_has_no_reads )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 2 ), "AA" );
    auto read1 = makeRead( startPos, startPos + 1 );

    auto first = std::make_shared< Variant >( reference, Region( "1", startPos, startPos + 1 ), "T" );

    auto second = std::make_shared< Variant >( reference, Region( "1", startPos + 1, startPos + 2 ), "T" );
    second->addRead( read1 );

    std::vector< readPtr_t > firstReads;
    std::vector< readPtr_t > secondReads;
    VariantCombinations( 1, 40 ).getOverlappingReads( first, second, firstReads, secondReads );
    BOOST_CHECK( firstReads.empty() );
    BOOST_CHECK( secondReads.empty() );
}

BOOST_AUTO_TEST_CASE( test_should_return_only_overlapping_reads_for_snps )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 3 ), "AAA" );

    auto read1 = makeRead( startPos, startPos + 2 );
    auto read2 = makeRead( startPos, startPos + 3 );
    auto read3 = makeRead( startPos + 1, startPos + 3 );

    auto first = std::make_shared< Variant >( reference, Region( "1", startPos, startPos + 1 ), "T" );
    first->addRead( read1 );
    first->addRead( read2 );

    auto second = std::make_shared< Variant >( reference, Region( "1", startPos + 2, startPos + 3 ), "T" );
    second->addRead( read2 );
    second->addRead( read3 );

    std::vector< readPtr_t > firstReads;
    std::vector< readPtr_t > secondReads;
    VariantCombinations( 1, 40 ).getOverlappingReads( first, second, firstReads, secondReads );

    BOOST_CHECK_EQUAL( firstReads.size(), 1 );
    BOOST_CHECK_EQUAL( firstReads[0], read2 );
    BOOST_CHECK_EQUAL( secondReads.size(), 1 );
    BOOST_CHECK_EQUAL( secondReads[0], read2 );
}

BOOST_AUTO_TEST_CASE( test_use_maximal_read_interval_in_reference )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 3 ), "AAA" );

    auto read1 = makeRead( startPos, startPos + 2 );
    auto read2 = std::make_shared< Read >( echidna::utils::BasePairSequence( 3, 'T' ), std::string( 3, 'Q' ), "",
                                           Cigar( "1S1M1S" ), 0, startPos + 1, 0, 0, 0, 0, 0, reference );
    auto read3 = std::make_shared< Read >( echidna::utils::BasePairSequence( 3, 'T' ), std::string( 3, 'Q' ), "",
                                           Cigar( "1I1M1I" ), 0, startPos + 1, 0, 0, 0, 0, 0, reference );
    auto read4 = makeRead( startPos + 1, startPos + 3 );

    auto first = std::make_shared< Variant >( reference, Region( "1", startPos, startPos + 1 ), "T" );
    first->addRead( read1 );
    first->addRead( read2 );
    first->addRead( read3 );

    auto second = std::make_shared< Variant >( reference, Region( "1", startPos + 2, startPos + 3 ), "T" );
    second->addRead( read2 );
    second->addRead( read3 );
    second->addRead( read4 );

    std::vector< readPtr_t > firstReads;
    std::vector< readPtr_t > secondReads;
    VariantCombinations( 1, 40 ).getOverlappingReads( first, second, firstReads, secondReads );

    BOOST_CHECK_EQUAL( firstReads.size(), 2 );
    BOOST_CHECK_EQUAL( firstReads[0], read2 );
    BOOST_CHECK_EQUAL( firstReads[1], read3 );
    BOOST_CHECK_EQUAL( secondReads.size(), 2 );
    BOOST_CHECK_EQUAL( secondReads[0], read2 );
    BOOST_CHECK_EQUAL( secondReads[1], read3 );
}

BOOST_AUTO_TEST_CASE( test_should_return_only_overlapping_reads_for_insertion_left_aligned_to_reference_start )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 5 ), "AAAAC" );

    auto read1 = makeRead( startPos + 1, startPos + 5 );  // overlaps only without left alignment
    auto read2 = makeRead( startPos, startPos + 5 );
    auto read3 = makeRead( startPos + 2, startPos + 5 );  // overlaps only without left alignment

    auto insertion = std::make_shared< Variant >( reference, Region( "1", startPos + 2, startPos + 2 ), "A" );
    insertion->addRead( read1 );
    insertion->addRead( read2 );
    insertion->addRead( read3 );

    auto snp = std::make_shared< Variant >( reference, Region( "1", startPos + 4, startPos + 5 ), "T" );
    snp->addRead( read1 );
    snp->addRead( read2 );
    snp->addRead( read3 );

    std::vector< readPtr_t > firstReads;
    std::vector< readPtr_t > secondReads;
    VariantCombinations( 1, 40 ).getOverlappingReads( insertion, snp, firstReads, secondReads );

    BOOST_CHECK_EQUAL( firstReads.size(), 1 );
    BOOST_CHECK_EQUAL( firstReads[0], read2 );
    BOOST_CHECK_EQUAL( secondReads.size(), 1 );
    BOOST_CHECK_EQUAL( secondReads[0], read2 );
}

BOOST_AUTO_TEST_CASE( test_should_return_only_overlapping_reads_for_insertion_left_aligned_to_reference )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 5 ), "TAAAC" );

    auto read1 = makeRead( startPos + 1, startPos + 5 );
    auto read2 = makeRead( startPos, startPos + 5 );
    auto read3 = makeRead( startPos + 2, startPos + 5 );  // overlaps only without left alignment

    auto insertion = std::make_shared< Variant >( reference, Region( "1", startPos + 2, startPos + 2 ), "A" );
    insertion->addRead( read1 );
    insertion->addRead( read2 );
    insertion->addRead( read3 );

    auto snp = std::make_shared< Variant >( reference, Region( "1", startPos + 4, startPos + 5 ), "T" );
    snp->addRead( read1 );
    snp->addRead( read2 );
    snp->addRead( read3 );

    std::vector< readPtr_t > firstReads;
    std::vector< readPtr_t > secondReads;
    VariantCombinations( 1, 40 ).getOverlappingReads( insertion, snp, firstReads, secondReads );

    BOOST_CHECK_EQUAL( firstReads.size(), 2 );
    BOOST_CHECK_EQUAL( firstReads[0], read1 );
    BOOST_CHECK_EQUAL( firstReads[1], read2 );
    BOOST_CHECK_EQUAL( secondReads.size(), 2 );
    BOOST_CHECK_EQUAL( secondReads[0], read1 );
    BOOST_CHECK_EQUAL( secondReads[1], read2 );
}

BOOST_AUTO_TEST_CASE( test_should_return_only_overlapping_reads_for_insertion_right_aligned_to_reference_end )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 5 ), "CAAAA" );

    auto read1 = makeRead( startPos, startPos + 3 );  // overlaps only without right alignment
    auto read2 = makeRead( startPos, startPos + 4 );  // overlaps only without right alignment
    auto read3 = makeRead( startPos, startPos + 5 );

    auto insertion = std::make_shared< Variant >( reference, Region( "1", startPos + 3, startPos + 3 ), "A" );
    insertion->addRead( read1 );
    insertion->addRead( read2 );
    insertion->addRead( read3 );

    auto snp = std::make_shared< Variant >( reference, Region( "1", startPos, startPos + 1 ), "T" );
    snp->addRead( read1 );
    snp->addRead( read2 );
    snp->addRead( read3 );

    std::vector< readPtr_t > firstReads;
    std::vector< readPtr_t > secondReads;
    VariantCombinations( 1, 40 ).getOverlappingReads( insertion, snp, firstReads, secondReads );

    BOOST_CHECK_EQUAL( firstReads.size(), 1 );
    BOOST_CHECK_EQUAL( firstReads[0], read3 );
    BOOST_CHECK_EQUAL( secondReads.size(), 1 );
    BOOST_CHECK_EQUAL( secondReads[0], read3 );
}

BOOST_AUTO_TEST_CASE( test_should_return_only_overlapping_reads_for_insertion_right_aligned_to_reference )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 5 ), "CAAAC" );

    auto read1 = makeRead( startPos, startPos + 3 );  // overlaps only without right alignment
    auto read2 = makeRead( startPos, startPos + 4 );
    auto read3 = makeRead( startPos, startPos + 5 );

    auto insertion = std::make_shared< Variant >( reference, Region( "1", startPos + 3, startPos + 3 ), "A" );
    insertion->addRead( read1 );
    insertion->addRead( read2 );
    insertion->addRead( read3 );

    auto snp = std::make_shared< Variant >( reference, Region( "1", startPos, startPos + 1 ), "T" );
    snp->addRead( read1 );
    snp->addRead( read2 );
    snp->addRead( read3 );

    std::vector< readPtr_t > firstReads;
    std::vector< readPtr_t > secondReads;
    VariantCombinations( 1, 40 ).getOverlappingReads( insertion, snp, firstReads, secondReads );

    BOOST_CHECK_EQUAL( firstReads.size(), 2 );
    BOOST_CHECK_EQUAL( firstReads[0], read2 );
    BOOST_CHECK_EQUAL( firstReads[1], read3 );
    BOOST_CHECK_EQUAL( secondReads.size(), 2 );
    BOOST_CHECK_EQUAL( secondReads[0], read2 );
    BOOST_CHECK_EQUAL( secondReads[1], read3 );
}

BOOST_AUTO_TEST_CASE( test_should_return_only_overlapping_reads_for_deletion_right_aligned_to_reference )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 4 ), "CACC" );

    auto read1 = makeRead( startPos, startPos + 2 );
    auto read2 = makeRead( startPos, startPos + 3 );  // overlaps only without right alignment
    auto read3 = makeRead( startPos, startPos + 4 );

    auto deletion = std::make_shared< Variant >( reference, Region( "1", startPos + 2, startPos + 3 ), "" );
    deletion->addRead( read1 );
    deletion->addRead( read2 );
    deletion->addRead( read3 );

    auto snp = std::make_shared< Variant >( reference, Region( "1", startPos, startPos + 1 ), "T" );
    snp->addRead( read1 );
    snp->addRead( read2 );
    snp->addRead( read3 );

    std::vector< readPtr_t > firstReads;
    std::vector< readPtr_t > secondReads;
    VariantCombinations( 1, 40 ).getOverlappingReads( deletion, snp, firstReads, secondReads );

    BOOST_CHECK_EQUAL( firstReads.size(), 1 );
    BOOST_CHECK_EQUAL( firstReads[0], read3 );
    BOOST_CHECK_EQUAL( secondReads.size(), 1 );
    BOOST_CHECK_EQUAL( secondReads[0], read3 );
}

BOOST_AUTO_TEST_CASE( test_should_return_only_overlapping_reads_for_deletion_left_aligned_to_reference )
{
    const int64_t startPos = 0;

    auto reference = std::make_shared< ReferenceSequence >( Region( "1", startPos, startPos + 4 ), "CCAC" );

    auto read1 = makeRead( startPos, startPos + 4 );
    auto read2 = makeRead( startPos + 1, startPos + 4 );  // overlaps only without right alignment
    auto read3 = makeRead( startPos + 2, startPos + 4 );

    auto deletion = std::make_shared< Variant >( reference, Region( "1", startPos + 1, startPos + 2 ), "" );
    deletion->addRead( read1 );
    deletion->addRead( read2 );
    deletion->addRead( read3 );

    auto snp = std::make_shared< Variant >( reference, Region( "1", startPos + 3, startPos + 4 ), "T" );
    snp->addRead( read1 );
    snp->addRead( read2 );
    snp->addRead( read3 );

    std::vector< readPtr_t > firstReads;
    std::vector< readPtr_t > secondReads;
    VariantCombinations( 1, 40 ).getOverlappingReads( deletion, snp, firstReads, secondReads );

    BOOST_CHECK_EQUAL( firstReads.size(), 1 );
    BOOST_CHECK_EQUAL( firstReads[0], read1 );
    BOOST_CHECK_EQUAL( secondReads.size(), 1 );
    BOOST_CHECK_EQUAL( secondReads[0], read1 );
}