#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "assembly/sequenceGraph.hpp"
#include "common.hpp"

using echidna::assembly::Node;
using echidna::assembly::SequenceGraph;
using echidna::assembly::SequenceGraph;
using echidna::utils::ReferenceSequence;
using echidna::caller::Region;
using echidna::utils::BasePairSequence;
using echidna::utils::QualitySequence;

BOOST_AUTO_TEST_CASE( testKmerSequenceConstruction )
{
    const BasePairSequence readSequence( "AKTCT" );

    const SequenceGraph::KmerSequence kmerSequence( &readSequence, 1, 4 );
    BOOST_CHECK_EQUAL( kmerSequence.getSequence(), readSequence.substr( 1, 4 ) );
    BOOST_CHECK_EQUAL( kmerSequence.front(), 'K' );
    BOOST_CHECK_EQUAL( kmerSequence.back(), 'T' );
}

BOOST_AUTO_TEST_CASE( testKmerSequenceSort )
{
    const BasePairSequence readSequence( "AATCT" );

    // Sorts on length first.
    BOOST_CHECK( SequenceGraph::KmerSequence( &readSequence, 0, 1 ) <
                 SequenceGraph::KmerSequence( &readSequence, 0, 2 ) );

    // Then sorts on sequence string.
    BOOST_CHECK( SequenceGraph::KmerSequence( &readSequence, 3, 1 ) <
                 SequenceGraph::KmerSequence( &readSequence, 2, 1 ) );
}

BOOST_AUTO_TEST_CASE( testKmerSequenceEquality )
{
    const BasePairSequence readSequence( "AATCT" );

    // Equality should fall to string equality.
    BOOST_CHECK( SequenceGraph::KmerSequence( &readSequence, 0, 1 ) ==
                 SequenceGraph::KmerSequence( &readSequence, 1, 1 ) );
}

BOOST_AUTO_TEST_CASE( testAddingReadSequenceWithFewerBPsThanKmerSize )
{
    SequenceGraph sequenceGraph( 10, 0 );
    const auto disallowRepeats = false;
    const BasePairSequence readSequence( "AATCA" );
    BOOST_CHECK_THROW( sequenceGraph.addReadSequence( &readSequence, QualitySequence( 5, 20 ), disallowRepeats ),
                       echidna::utils::echidna_exception );
}

BOOST_AUTO_TEST_CASE( testAddingReadSequenceContructsTwoNodesWithEdge )
{
    SequenceGraph sequenceGraph( 4, 0 );
    const BasePairSequence readSequence( "AATCA" );
    sequenceGraph.addReadSequence( &readSequence, QualitySequence( 5, 20 ), 0 );

    BOOST_CHECK( sequenceGraph.containsNode( "AATC" ) );
    BOOST_CHECK( sequenceGraph.containsNode( "ATCA" ) );
}

BOOST_AUTO_TEST_CASE( testAddingReferenceSequenceContructsTwoNodesWithEdge )
{
    SequenceGraph sequenceGraph( 4, 0 );

    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "AATCA" );
    sequenceGraph.addReferenceSequence( referenceSequence, referenceSequence->region(), false );

    BOOST_CHECK( sequenceGraph.containsNode( "AATC" ) );
    BOOST_CHECK( sequenceGraph.containsNode( "ATCA" ) );
}

BOOST_AUTO_TEST_CASE( testAddingReferenceShouldReturnWhetherARepeatedPositionHasBeenAdded )
{
    const auto referenceSequence0 = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "AATCA" );
    const auto referenceSequence1 = std::make_shared< ReferenceSequence >( Region( "1", 0, 8 ), "AATCAATC" );

    BOOST_CHECK(
        not SequenceGraph( 4, 0 ).addReferenceSequence( referenceSequence0, referenceSequence0->region(), false ) );
    BOOST_CHECK(
        SequenceGraph( 4, 0 ).addReferenceSequence( referenceSequence1, referenceSequence1->region(), false ) );
}

BOOST_AUTO_TEST_CASE( testShouldNotAddNodeIfRepeatsAreDisallowedAndReturnEarly )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 9 ), "AATCAATCG" );
    const std::size_t kmerSize = 4;

    {
        SequenceGraph sequenceGraph( kmerSize, 0 );
        bool disallowRepeats = false;
        BOOST_CHECK(
            sequenceGraph.addReferenceSequence( referenceSequence, referenceSequence->region(), disallowRepeats ) );
        BOOST_CHECK_EQUAL( sequenceGraph.size(), 5 );
    }
    {
        SequenceGraph sequenceGraph( kmerSize, 0 );
        bool disallowRepeats = true;
        BOOST_CHECK(
            sequenceGraph.addReferenceSequence( referenceSequence, referenceSequence->region(), disallowRepeats ) );
        BOOST_CHECK_EQUAL( sequenceGraph.size(), 4 );
    }
}

BOOST_AUTO_TEST_CASE( testShouldNotAddNodeIfRepeatsAreDisallowedAndReturnEarlyForReadSequence )
{
    const auto readSequence = BasePairSequence( "AATCAATCG" );
    const auto qualitySeque = "QQQQQQQQQ";
    const std::size_t kmerSize = 4;

    {
        SequenceGraph sequenceGraph( kmerSize, 0 );
        bool disallowRepeats = false;
        BOOST_CHECK( sequenceGraph.addReadSequence( &readSequence, qualitySeque, disallowRepeats ) );
        BOOST_CHECK_EQUAL( sequenceGraph.size(), 5 );
    }
    {
        SequenceGraph sequenceGraph( kmerSize, 0 );
        bool disallowRepeats = true;
        BOOST_CHECK( sequenceGraph.addReadSequence( &readSequence, qualitySeque, disallowRepeats ) );
        BOOST_CHECK_EQUAL( sequenceGraph.size(), 4 );
    }
}

BOOST_AUTO_TEST_CASE( testAddingReadSequenceShouldReturnWhetherARepeatedNodeHasBeenAdded )
{
    const auto disallowRepeats = false;
    const auto readSequence0 = BasePairSequence( "AATCA" );
    const auto readSequence1 = BasePairSequence( "AATCAATC" );
    BOOST_CHECK(
        not SequenceGraph( 4, 0 ).addReadSequence( &readSequence0, QualitySequence( 5, 20 ), disallowRepeats ) );
    BOOST_CHECK( SequenceGraph( 4, 0 ).addReadSequence( &readSequence1, QualitySequence( 8, 20 ), disallowRepeats ) );
}

BOOST_AUTO_TEST_CASE( testShouldGetBranchPointWithMultipleInNodes )
{
    SequenceGraph sequenceGraph( 4, 0 );
    const auto disallowRepeats = false;
    const auto readSequence0 = BasePairSequence( "AATTTTC" );
    const auto readSequence1 = BasePairSequence( "CATTTTC" );
    sequenceGraph.addReadSequence( &readSequence0, QualitySequence( 7, 20 ), disallowRepeats );
    sequenceGraph.addReadSequence( &readSequence1, QualitySequence( 7, 20 ), disallowRepeats );

    const auto chains = sequenceGraph.getChains( 0 );
    BOOST_REQUIRE_EQUAL( chains.size(), 3 );
    for ( const auto & chain : chains )
    {
        BOOST_CHECK( chain.front()->isBranch() );
        BOOST_CHECK( chain.back()->isBranch() );
    }
}

BOOST_AUTO_TEST_CASE( testChainSupportIsMinQualityAlongChain )
{
    const std::size_t minBaseQuality = 20;
    SequenceGraph sequenceGraph( 4, minBaseQuality );
    const auto disallowRepeats = false;
    const auto readSequence0 = BasePairSequence( "AATTTTC" );
    const auto readSequence1 = BasePairSequence( "CATTTTC" );
    const auto readSequence2 = BasePairSequence( "CATTTTC" );
    sequenceGraph.addReadSequence( &readSequence0, QualitySequence( 7, 20 ), disallowRepeats );
    sequenceGraph.addReadSequence( &readSequence1, QualitySequence( 7, 20 ), disallowRepeats );
    sequenceGraph.addReadSequence( &readSequence2, QualitySequence( 7, 19 ), disallowRepeats );

    const auto chains = sequenceGraph.getChains( 0 );

    BOOST_REQUIRE_EQUAL( chains.size(), 3 );

    BOOST_CHECK_EQUAL( chains[0].initialSequence(), "AATT" );
    BOOST_CHECK_EQUAL( chains[0].endSequence(), "ATTT" );
    BOOST_CHECK_EQUAL( chains[0].getSupport(), 20 );

    BOOST_CHECK_EQUAL( chains[1].initialSequence(), "ATTT" );
    BOOST_CHECK_EQUAL( chains[1].endSequence(), "TTTC" );
    BOOST_CHECK_EQUAL( chains[1].getSupport(), 40 );

    BOOST_CHECK_EQUAL( chains[2].initialSequence(), "CATT" );
    BOOST_CHECK_EQUAL( chains[2].endSequence(), "ATTT" );
    BOOST_CHECK_EQUAL( chains[2].getSupport(), 20 );
}

BOOST_AUTO_TEST_CASE( testShouldGetChainForInnerLoop )
{
    SequenceGraph sequenceGraph( 4, 0 );
    const auto disallowRepeats = false;
    const auto readSequence0 = BasePairSequence( "ATTTTTC" );
    sequenceGraph.addReadSequence( &readSequence0, QualitySequence( 7, 20 ), disallowRepeats );
    const auto chains = sequenceGraph.getChains( 0 );

    BOOST_REQUIRE_EQUAL( chains.size(), 3 );

    for ( const auto & chain : chains )
    {
        BOOST_CHECK( chain.front()->isBranch() );
        BOOST_CHECK( chain.back()->isBranch() );
        BOOST_CHECK_EQUAL( chain.size(), 2 );
    }

    BOOST_CHECK_EQUAL( chains[0].initialSequence(), "ATTT" );
    BOOST_CHECK_EQUAL( chains[0].endSequence(), "TTTT" );

    BOOST_CHECK_EQUAL( chains[1].initialSequence(), "TTTT" );
    BOOST_CHECK_EQUAL( chains[1].endSequence(), "TTTC" );

    // Inner loop in graph starts and ends at same node.
    BOOST_CHECK_EQUAL( chains[2].initialSequence(), "TTTT" );
    BOOST_CHECK_EQUAL( chains[2].endSequence(), "TTTT" );
}

BOOST_AUTO_TEST_CASE( testShouldRemoveCycle )
{
    SequenceGraph sequenceGraph( 2, 0 );
    const auto disallowRepeats = false;
    const auto readSequence0 = BasePairSequence( "AATTAA" );
    sequenceGraph.addReadSequence( &readSequence0, QualitySequence( 6, 20 ), disallowRepeats );

    // Should a pure cycle be classified as a chain? This would mean all edges belong to exactly one chain.
    // The only problem is that the start end is not well-defined but with a chain class this wouldn't be an issue.
    //    const auto chains = sequenceGraph.getChains();
    //    BOOST_REQUIRE_EQUAL( chains.size(), 1 );

    // Once the first thing is checked should also check that the cycle can be removed.
}
