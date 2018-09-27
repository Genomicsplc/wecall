// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "variant/type/variant.hpp"

#include "variant/clustering.hpp"

using echidna::variant::variantSet_t;
using echidna::utils::ReferenceSequence;
using echidna::caller::Region;
using echidna::variant::VariantCluster;
using echidna::variant::Variant;
using echidna::variant::varPtr_t;

using namespace echidna::alignment;
using namespace echidna::io;
using echidna::utils::ReferenceSequence;

BOOST_AUTO_TEST_CASE( testComputationOfPaddingRegionsForEmptyCluster )
{
    const Region blockRegion( "1", 0, 100 );
    BOOST_CHECK_EQUAL( echidna::variant::computeClustersPaddingRegions( blockRegion, {} ).size(), 0 );
}

BOOST_AUTO_TEST_CASE( testComputationOfPaddingRegionsForOneCluster )
{
    const Region blockRegion( "1", 0, 100 );

    std::vector< VariantCluster > clusters;
    clusters.emplace_back( std::vector< varPtr_t >{}, Region( "1", 1, 99 ) );

    const auto regions = echidna::variant::computeClustersPaddingRegions( blockRegion, clusters );
    BOOST_CHECK_EQUAL( regions.size(), 1 );
    BOOST_CHECK_EQUAL( regions.front(), Region( "1", 0, 100 ) );
}

BOOST_AUTO_TEST_CASE( testComputationOfPaddingRegionsForManyClusters )
{
    const Region blockRegion( "1", 0, 100 );

    std::vector< VariantCluster > clusters;
    clusters.emplace_back( std::vector< varPtr_t >{}, Region( "1", 10, 40 ) );
    clusters.emplace_back( std::vector< varPtr_t >{}, Region( "1", 60, 90 ) );

    const auto regions = echidna::variant::computeClustersPaddingRegions( blockRegion, clusters );
    BOOST_REQUIRE_EQUAL( regions.size(), 2 );
    BOOST_CHECK_EQUAL( regions[0], Region( "1", 0, 60 ) );
    BOOST_CHECK_EQUAL( regions[1], Region( "1", 40, 100 ) );
}

BOOST_AUTO_TEST_CASE( testStartEndOfVariantClusterSNPAndMNP )
{  // 01234
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 15 ), "AAAAA" );
    auto mnp = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "TTTTT" );
    auto snp = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "T" );
    std::vector< echidna::variant::varPtr_t > variants = {mnp, snp};

    std::sort( variants.begin(), variants.end(), echidna::variant::varPtrComp() );

    VariantCluster cluster( variants, echidna::caller::Region( "1", 10, 15 ) );

    BOOST_CHECK_EQUAL_COLLECTIONS( cluster.variants().begin(), cluster.variants().end(), variants.begin(),
                                   variants.end() );
    BOOST_CHECK_EQUAL( cluster.region().start(), 10 );
    BOOST_CHECK_EQUAL( cluster.zeroIndexedVcfStart( {} ), 10 );
    BOOST_CHECK_EQUAL( cluster.region().end(), 15 );
}

BOOST_AUTO_TEST_CASE( testStartEndOfVariantClusterSNPAndDel )
{  // 012345
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 11, 16 ), "AAAAA" );
    auto del = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 16 ), "" );
    auto snp = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "T" );
    std::vector< echidna::variant::varPtr_t > variants = {del, snp};

    std::sort( variants.begin(), variants.end(), echidna::variant::varPtrComp() );

    VariantCluster cluster( variants, echidna::caller::Region( "1", 11, 16 ) );

    echidna::caller::Call call( del, del->interval(), 100, 1, {{}} );

    BOOST_CHECK_EQUAL_COLLECTIONS( cluster.variants().begin(), cluster.variants().end(), variants.begin(),
                                   variants.end() );
    BOOST_CHECK_EQUAL( cluster.region().start(), 11 );
    BOOST_CHECK_EQUAL( cluster.zeroIndexedVcfStart( {call} ), 10 );
    BOOST_CHECK_EQUAL( cluster.region().end(), 16 );
}

BOOST_AUTO_TEST_CASE( testStartEndOfVariantClusterDelAndMNP )
{  // 012345
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 13 ), "AAA" );
    auto del = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 13 ), "" );
    auto mnp = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "TTT" );
    std::vector< echidna::variant::varPtr_t > variants = {del, mnp};

    std::sort( variants.begin(), variants.end(), echidna::variant::varPtrComp() );

    VariantCluster cluster( variants, echidna::caller::Region( "1", 10, 13 ) );

    BOOST_CHECK_EQUAL_COLLECTIONS( cluster.variants().begin(), cluster.variants().end(), variants.begin(),
                                   variants.end() );
    BOOST_CHECK_EQUAL( cluster.region().start(), 10 );
    BOOST_CHECK_EQUAL( cluster.zeroIndexedVcfStart( {} ), 10 );
    BOOST_CHECK_EQUAL( cluster.region().end(), 13 );
}

BOOST_AUTO_TEST_CASE( testVariantClusterGenerationSNPandMNP )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 22 ), "AAAAAAAAAAAA" );
    auto snp1 = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 11 ), "C" );
    // 012345678901
    auto longMnp = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "TTTTTTTAATTC" );
    auto snp2 = std::make_shared< Variant >( referenceSequence, Region( "1", 21, 22 ), "T" );

    variantSet_t variants = {snp1, longMnp, snp2};

    int minClusterDistance = 0;

    auto clusters = echidna::variant::generateVariantClusters( variants, minClusterDistance, Region( "1", 10, 25 ) );
    BOOST_CHECK_EQUAL( clusters.size(), 1 );
    BOOST_CHECK_EQUAL_COLLECTIONS( clusters[0].variants().begin(), clusters[0].variants().end(), variants.begin(),
                                   variants.end() );
}

BOOST_AUTO_TEST_CASE( testVariantClusterGenerationSNPandDel )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 22 ), "AAAAAAAAAAAA" );
    // 012345
    auto del = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 16 ), "" );
    auto snp = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "T" );
    auto snp2 = std::make_shared< Variant >( referenceSequence, Region( "1", 17, 18 ), "T" );

    variantSet_t variants = {del, snp, snp2};

    int minClusterDistance = 0;

    auto clusters = echidna::variant::generateVariantClusters( variants, minClusterDistance, Region( "1", 10, 20 ) );
    BOOST_CHECK_EQUAL( clusters.size(), 1 );
    BOOST_CHECK_EQUAL_COLLECTIONS( clusters[0].variants().begin(), clusters[0].variants().end(), variants.begin(),
                                   variants.end() );
}

BOOST_AUTO_TEST_CASE( testVariantClusterGenerationSNPAndMNP )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 5, 15 ), "AAAAAAAAAA" );
    auto snp = std::make_shared< Variant >( referenceSequence, Region( "1", 6, 7 ), "T" );
    // 012345
    auto mnp = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 13 ), "TTT" );
    variantSet_t variants = {snp, mnp};

    int minClusterDistance = 1;
    int maxClusterDistance = 4;
    int maxClusterSize = 200;

    const auto clusters =
        echidna::variant::generateMergedClusters( variants, minClusterDistance, maxClusterDistance, maxClusterSize,
                                                  referenceSequence->region(), referenceSequence, 100, 0 );
    BOOST_CHECK_EQUAL( clusters.size(), 1 );
    BOOST_CHECK_EQUAL_COLLECTIONS( clusters[0].variants().begin(), clusters[0].variants().end(), variants.begin(),
                                   variants.end() );
}

BOOST_AUTO_TEST_CASE( testVariantClusterGenerationSNPAndDelAndMNP )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 5, 15 ), "AAAAAAAAAA" );
    auto snp = std::make_shared< Variant >( referenceSequence, Region( "1", 6, 7 ), "T" );
    // 012345
    auto del = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 13 ), "" );
    auto mnp = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 13 ), "TTT" );
    variantSet_t variants = {snp, del, mnp};

    int minClusterDistance = 0;

    auto clusters =
        echidna::variant::generateVariantClusters( variants, minClusterDistance, referenceSequence->region() );

    BOOST_CHECK_EQUAL( clusters.size(), 1 );
    BOOST_CHECK_EQUAL_COLLECTIONS( clusters[0].variants().begin(), clusters[0].variants().end(), variants.begin(),
                                   variants.end() );
}

BOOST_AUTO_TEST_CASE( testVariantsOverlappingAnAlignmentOfIndelGetAddedToSameCluster )
{
    auto referenceSequence =
        std::make_shared< ReferenceSequence >( Region( "1", 0, 100 ), echidna::utils::BasePairSequence( 100, 'A' ) );

    auto snp1 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "T" );
    auto snp2 = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "T" );

    auto snp3 = std::make_shared< Variant >( referenceSequence, Region( "1", 99, 100 ), "T" );

    auto del = std::make_shared< Variant >( referenceSequence, Region( "1", 50, 51 ), "" );

    variantSet_t variants = {snp1, snp2, del, snp3};

    int minClusterDistance = 0;

    const auto clusters =
        echidna::variant::generateVariantClusters( variants, minClusterDistance, referenceSequence->region() );

    BOOST_CHECK_EQUAL( clusters.size(), 1 );

    BOOST_CHECK_EQUAL_COLLECTIONS( clusters[0].variants().begin(), clusters[0].variants().end(), variants.begin(),
                                   variants.end() );
}

BOOST_AUTO_TEST_CASE( testShouldNotMergeClustersIfVariantCombinationsNotComputed )
{
    const auto reference =
        std::make_shared< ReferenceSequence >( ReferenceSequence( Region( "1", 1, 12 ), "AAAAAAAAAAA" ) );
    VariantCluster cluster1( {std::make_shared< Variant >( reference, Region( "1", 1, 2 ), "T" )},
                             Region( "1", 1, 2 ) );
    VariantCluster cluster2( {std::make_shared< Variant >( reference, Region( "1", 11, 12 ), "T" )},
                             Region( "1", 11, 12 ) );

    const auto maxCombinations = 4;
    const auto mergedClusters = echidna::variant::mergeClusters( {cluster1, cluster2}, maxCombinations, 9, 100 );

    BOOST_CHECK_EQUAL( mergedClusters.size(), 2 );
}

BOOST_AUTO_TEST_CASE( testMergeClusters )
{
    const auto reference =
        std::make_shared< ReferenceSequence >( ReferenceSequence( Region( "1", 1, 12 ), "AAAAAAAAAAA" ) );
    VariantCluster cluster1( {std::make_shared< Variant >( reference, Region( "1", 1, 2 ), "T" )},
                             Region( "1", 1, 2 ) );
    VariantCluster cluster2( {std::make_shared< Variant >( reference, Region( "1", 11, 12 ), "T" )},
                             Region( "1", 11, 12 ) );

    const auto maxCombinations = 4;
    cluster1.computeVariantCombinations( 5, 40, maxCombinations, reference );
    cluster2.computeVariantCombinations( 5, 40, maxCombinations, reference );
    const auto mergedClusters = echidna::variant::mergeClusters( {cluster1, cluster2}, maxCombinations, 9, 100 );

    BOOST_CHECK_EQUAL( mergedClusters.size(), 1 );

    BOOST_CHECK_EQUAL( mergedClusters[0].variants().size(), 2 );
}

BOOST_AUTO_TEST_CASE( testMaxVariantCombinationsGetsTwoClusters )
{
    const auto reference =
        std::make_shared< ReferenceSequence >( ReferenceSequence( Region( "1", 1, 12 ), "AAAAAAAAAAA" ) );
    VariantCluster cluster1( {std::make_shared< Variant >( reference, Region( "1", 1, 2 ), "T" )},
                             Region( "1", 1, 2 ) );
    VariantCluster cluster2( {std::make_shared< Variant >( reference, Region( "1", 11, 12 ), "T" )},
                             Region( "1", 11, 12 ) );

    const auto maxCombinations = 2;
    cluster1.computeVariantCombinations( 1, 40, maxCombinations, reference );
    cluster2.computeVariantCombinations( 1, 40, maxCombinations, reference );
    const auto mergedClusters = echidna::variant::mergeClusters( {cluster1, cluster2}, maxCombinations, 9, 100 );

    BOOST_CHECK_EQUAL( mergedClusters.size(), 2 );

    BOOST_CHECK_EQUAL( mergedClusters[0].variants().size(), 1 );
    BOOST_CHECK_EQUAL( mergedClusters[1].variants().size(), 1 );
}

BOOST_AUTO_TEST_CASE( testMaxVariantsDistance )
{
    const auto reference =
        std::make_shared< ReferenceSequence >( ReferenceSequence( Region( "1", 1, 12 ), "AAAAAAAAAAA" ) );
    VariantCluster cluster1( {std::make_shared< Variant >( reference, Region( "1", 1, 2 ), "T" )},
                             Region( "1", 1, 2 ) );
    VariantCluster cluster2( {std::make_shared< Variant >( reference, Region( "1", 11, 12 ), "T" )},
                             Region( "1", 11, 12 ) );

    const auto maxCombinations = 5;
    cluster1.computeVariantCombinations( 2, 40, maxCombinations, reference );
    cluster2.computeVariantCombinations( 2, 40, maxCombinations, reference );
    const auto mergedClusters = echidna::variant::mergeClusters( {cluster1, cluster2}, maxCombinations, 8, 100 );

    BOOST_CHECK_EQUAL( mergedClusters.size(), 2 );

    BOOST_CHECK_EQUAL( mergedClusters[0].variants().size(), 1 );
    BOOST_CHECK_EQUAL( mergedClusters[1].variants().size(), 1 );
}

BOOST_AUTO_TEST_CASE( testCreationOfSubclusters )
{
    const auto reference =
        std::make_shared< ReferenceSequence >( ReferenceSequence( Region( "1", 99, 120 ), "ATAAAAAAAAAAAAAAAAAAC" ) );
    VariantCluster cluster1( {std::make_shared< Variant >( reference, Region( "1", 99, 100 ), "C" ),
                              std::make_shared< Variant >( reference, Region( "1", 100, 120 ), "T" ),
                              std::make_shared< Variant >( reference, Region( "1", 105, 105 ), "T" ),
                              std::make_shared< Variant >( reference, Region( "1", 110, 110 ), "T" ),
                              std::make_shared< Variant >( reference, Region( "1", 115, 120 ), "T" )},
                             reference->region() );

    const auto subclusters = cluster1.buildSubClusters( 4 );

    const auto mainCluster = std::get< 0 >( subclusters );
    const auto miniClusters = std::get< 1 >( subclusters );

    BOOST_CHECK_EQUAL( mainCluster.region(), reference->region() );
    BOOST_REQUIRE_EQUAL( mainCluster.variants().size(), 3 );
    BOOST_CHECK_EQUAL( *mainCluster.variants()[0], *cluster1.variants()[0] );
    BOOST_CHECK_EQUAL( *mainCluster.variants()[1], *cluster1.variants()[1] );
    BOOST_CHECK_EQUAL( *mainCluster.variants()[2], *cluster1.variants()[4] );
    BOOST_CHECK_EQUAL( mainCluster.region(), Region( "1", 99, 120 ) );

    BOOST_REQUIRE_EQUAL( miniClusters.size(), 2 );

    BOOST_REQUIRE_EQUAL( miniClusters[0].variants().size(), 1 );
    BOOST_CHECK_EQUAL( *miniClusters[0].variants()[0], *cluster1.variants()[2] );
    BOOST_CHECK_EQUAL( miniClusters[0].region(), Region( "1", 105, 105 ) );
    BOOST_CHECK_EQUAL( miniClusters[0].paddedRegion(), Region( "1", 100, 110 ) );

    BOOST_REQUIRE_EQUAL( miniClusters[1].variants().size(), 1 );
    BOOST_CHECK_EQUAL( *miniClusters[1].variants()[0], *cluster1.variants()[3] );
    BOOST_CHECK_EQUAL( miniClusters[1].region(), Region( "1", 110, 110 ) );
    BOOST_CHECK_EQUAL( miniClusters[1].paddedRegion(), Region( "1", 105, 115 ) );
}
