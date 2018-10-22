// All content Copyright (C) 2018 Genomics plc
#include "variant/type/variant.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <string>

#include "utils/exceptions.hpp"
#include "utils/referenceSequence.hpp"

using wecall::caller::Region;
using wecall::caller::SetRegions;
using wecall::utils::ReferenceSequence;
using wecall::variant::Variant;

BOOST_AUTO_TEST_CASE( testInsertion )
{
    using namespace wecall::variant;

    const std::string contig( "1" );
    const std::string added( "A" );
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( contig, 1000001, 1000001 ), "" );
    varPtr_t theInsertion( new Variant( referenceSequence, referenceSequence->region(), added, false ) );

    BOOST_CHECK_EQUAL( theInsertion->sequence(), std::string( "A" ) );
    BOOST_CHECK_EQUAL( theInsertion->start(), 1000001 );
    BOOST_CHECK_EQUAL( theInsertion->end(), 1000001 );

    auto expectedPrior = 1e-4 * pow( 0.33, added.size() );
    wecall::variant::setDefaultPriors( {theInsertion} );
    BOOST_CHECK_CLOSE( theInsertion->prior(), expectedPrior, 1e-5 );
}

BOOST_AUTO_TEST_CASE( testLeftAlignIns )
{
    using namespace wecall::variant;

    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 999999, 1000003 ), "AGGG" );

    // const std::string added("G");
    varPtr_t theInsertion( new Variant( refSeq, Region( "1", 1000003, 1000003 ), "G", false ) );

    // Check that theInsertion left-aligns correctly
    theInsertion = theInsertion->getLeftAligned( 999999 );
    BOOST_CHECK_EQUAL( theInsertion->sequence(), std::string( "G" ) );
    BOOST_CHECK_EQUAL( theInsertion->start(), 1000000 );
    BOOST_CHECK_EQUAL( theInsertion->end(), 1000000 );
}

BOOST_AUTO_TEST_CASE( testLeftAlignIns2 )
{
    // Test left-aligning insertions
    using namespace wecall::variant;

    // Reference genome 37, chrom 1 from 1000000 - 1000099
    //
    //                          0         1         2         3         4         5         6         7         8 9
    //                          0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    const auto refSeq = std::make_shared< ReferenceSequence >(
        Region( "1", 1000000, 1000100 ),
        "GGGCACAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACGGGAAGGAA" );

    varPtr_t ins1( new Variant( refSeq, Region( "1", 1000054, 1000054 ), "C", false ) );
    varPtr_t ins2( new Variant( refSeq, Region( "1", 1000058, 1000058 ), "C", false ) );
    varPtr_t ins3( new Variant( refSeq, Region( "1", 1000058, 1000058 ), "CC", false ) );
    varPtr_t ins4( new Variant( refSeq, Region( "1", 1000058, 1000058 ), "CCC", false ) );
    varPtr_t ins5( new Variant( refSeq, Region( "1", 1000058, 1000058 ), "CCCCCCCCCCCC", false ) );
    varPtr_t ins6( new Variant( refSeq, Region( "1", 1000058, 1000058 ), "CCCATGCATGC", false ) );
    varPtr_t ins7( new Variant( refSeq, Region( "1", 1000055, 1000055 ), "A", false ) );

    // Chec ins1 left aligns correctly
    ins1 = ins1->getLeftAligned( 1000000 );
    BOOST_CHECK_EQUAL( ins1->sequence(), std::string( "C" ) );
    BOOST_CHECK_EQUAL( ins1->start(), 1000053 );
    BOOST_CHECK_EQUAL( ins1->end(), 1000053 );

    // Chec ins2 left aligns correctly
    ins2 = ins2->getLeftAligned( 1000000 );
    BOOST_CHECK_EQUAL( ins2->sequence(), std::string( "C" ) );
    BOOST_CHECK_EQUAL( ins2->start(), 1000053 );
    BOOST_CHECK_EQUAL( ins2->end(), 1000053 );

    // Chec ins3 left aligns correctly
    ins3 = ins3->getLeftAligned( 1000000 );
    BOOST_CHECK_EQUAL( ins3->sequence(), std::string( "CC" ) );
    BOOST_CHECK_EQUAL( ins3->start(), 1000053 );
    BOOST_CHECK_EQUAL( ins3->end(), 1000053 );

    // Chec ins4 left aligns correctly
    ins4 = ins4->getLeftAligned( 1000000 );
    BOOST_CHECK_EQUAL( ins4->sequence(), std::string( "CCC" ) );
    BOOST_CHECK_EQUAL( ins4->start(), 1000053 );
    BOOST_CHECK_EQUAL( ins4->end(), 1000053 );

    // Chec ins5 left aligns correctly
    ins5 = ins5->getLeftAligned( 1000000 );
    BOOST_CHECK_EQUAL( ins5->sequence(), std::string( "CCCCCCCCCCCC" ) );
    BOOST_CHECK_EQUAL( ins5->start(), 1000053 );
    BOOST_CHECK_EQUAL( ins5->end(), 1000053 );

    // Chec ins6 left aligns correctly
    ins6 = ins6->getLeftAligned( 1000000 );
    BOOST_CHECK_EQUAL( ins6->sequence(), std::string( "CCCCATGCATG" ) );
    BOOST_CHECK_EQUAL( ins6->start(), 1000057 );
    BOOST_CHECK_EQUAL( ins6->end(), 1000057 );

    // Chec ins7 left aligns correctly
    ins7 = ins7->getLeftAligned( 1000000 );
    BOOST_CHECK_EQUAL( ins7->sequence(), std::string( "A" ) );
    BOOST_CHECK_EQUAL( ins7->start(), 1000055 );
    BOOST_CHECK_EQUAL( ins7->end(), 1000055 );
}

BOOST_AUTO_TEST_CASE( insertionWontLeftAlignWhenBoundaryInsideRepetitiveRegion )
{
    // Test left-aligning insertions
    using namespace wecall::variant;

    // Reference genome 37, chrom 1 from 1000000 - 1000099
    //
    //                          0         1         2         3         4         5         6         7         8 9
    //                          0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    const auto refSeq = std::make_shared< ReferenceSequence >(
        Region( "1", 1000000, 1000100 ),
        "GGGCACAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACGGGAAGGAA" );

    varPtr_t ins( new Variant( refSeq, Region( "1", 1000057, 1000057 ), "C", false ) );

    // Chec ins1 left aligns correctly
    ins = ins->getLeftAligned( 1000053 );
    BOOST_CHECK_EQUAL( ins->start(), 1000053 );
    BOOST_CHECK_EQUAL( ins->end(), 1000053 );
}

BOOST_AUTO_TEST_CASE( insertionWillLeftAlignWhenBoundaryInsideRepetitiveRegion )
{
    // Test left-aligning insertions
    using namespace wecall::variant;

    // Reference genome 37, chrom 1 from 1000000 - 1000099
    //
    //                          0         1         2         3         4         5         6         7         8 9
    //                          0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    const auto refSeq = std::make_shared< ReferenceSequence >(
        Region( "1", 1000000, 1000100 ),
        "GGGCACAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACGGGAAGGAA" );

    varPtr_t ins( new Variant( refSeq, Region( "1", 1000057, 1000057 ), "C", false ) );

    // Chec ins1 left aligns correctly
    ins = ins->getLeftAligned( 1000052 );
    BOOST_CHECK_EQUAL( ins->start(), 1000053 );
}

BOOST_AUTO_TEST_CASE( insertionWillLeftAlignInComplexRepeatRegion )
{
    using namespace wecall::variant;
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 51, 101 ),
                                                               "TGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA" );
    varPtr_t ins( new Variant( refSeq, Region( "1", 101, 101 ), "GATTACA", false ) );

    ins = ins->getLeftAligned( 51 );
    BOOST_CHECK_EQUAL( ins->start(), 52 );
    BOOST_CHECK_EQUAL( ins->sequence(), std::string( "GATTACA" ) );
}

BOOST_AUTO_TEST_CASE( insertionWillLeftAlignInComplexRepeatRegionWithAllRefStringOfRepeatUnits )
{
    using namespace wecall::variant;
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 52, 101 ),
                                                               "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA" );
    varPtr_t ins( new Variant( refSeq, Region( "1", 101, 101 ), "GATTACA", false ) );

    ins = ins->getLeftAligned( 52 );
    BOOST_CHECK_EQUAL( ins->start(), 52 );
    BOOST_CHECK_EQUAL( ins->sequence(), std::string( "GATTACA" ) );
    BOOST_CHECK( not ins->isFullyLeftAligned() );
}

BOOST_AUTO_TEST_CASE( insertionWillLeftAlignInComplexRepeatRegionWithChangeToSequence )
{
    using namespace wecall::variant;
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( "1", 50, 101 ), "TAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA" );
    varPtr_t var( new Variant( referenceSequence, Region( "1", 101, 101 ), "GATTACA", false ) );
    //                                     0      1      2      3      4      5      6

    auto var2 = var->getLeftAligned( 50 );
    BOOST_CHECK_EQUAL( var2->start(), 51 );
    BOOST_CHECK_EQUAL( var2->sequence(), std::string( "AGATTAC" ) );
    BOOST_CHECK( var2->isFullyLeftAligned() );

    const auto var3 = var2->getRightAligned( var->end() );

    BOOST_CHECK_EQUAL( var->refSequence(), var3->refSequence() );
    BOOST_CHECK_EQUAL( var->sequence(), var3->sequence() );
    wecall::variant::setDefaultPriors( {var, var3} );
    BOOST_CHECK_CLOSE( var->prior(), var3->prior(), 1e-5 );
    BOOST_CHECK( not var3->isFullyLeftAligned() );
}

BOOST_AUTO_TEST_CASE(
    insertionWillLeftAlignInComplexRepeatRegionWithChangeToSequenceAndAllRefStringMatchingRepeatUnits )
{
    using namespace wecall::variant;
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( "1", 50, 101 ), "CAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA" );
    varPtr_t var = std::make_shared< Variant >( referenceSequence, Region( "1", 101, 101 ), "GATTACA", false );
    //                                     0      1      2      3      4      5      6

    auto var2 = var->getLeftAligned( 50 );
    BOOST_CHECK_EQUAL( var2->start(), 50 );
    BOOST_CHECK_EQUAL( var2->sequence(), std::string( "CAGATTA" ) );

    const auto var3 = var2->getRightAligned( var->end() );

    BOOST_CHECK_EQUAL( var->refSequence(), var3->refSequence() );
    BOOST_CHECK_EQUAL( var->sequence(), var3->sequence() );
    wecall::variant::setDefaultPriors( {var, var3} );
    BOOST_CHECK_CLOSE( var->prior(), var3->prior(), 1e-5 );
}

BOOST_AUTO_TEST_CASE( testInsertionJoinsWithDeletionOnlyOnRight )
{
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 10, 20 ), "CAGATTACAG" );

    const auto ins1 = std::make_shared< Variant >( refSeq, Region( "1", 10, 10 ), "A" );
    const auto ins2 = std::make_shared< Variant >( refSeq, Region( "1", 10, 10 ), "T" );
    const auto ins3 = std::make_shared< Variant >( refSeq, Region( "1", 11, 11 ), "C" );

    BOOST_REQUIRE( ins1->joinable( ins2 ) );

    const auto join = ins1->join( ins2 );

    BOOST_CHECK_EQUAL( join->region(), Region( "1", 10, 10 ) );
    BOOST_CHECK_EQUAL( join->sequence(), "AT" );

    BOOST_CHECK( not ins1->joinable( ins3 ) );
    BOOST_CHECK_THROW( ins1->join( ins3 ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( testInsertionBah )
{
    // 20:56258600-56258650 20:56258612-56258612
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "20", 56258600, 56258650 ),
                                                               "GCCTGTGTGCGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTTTCTTTTAAA" );
    const auto ins = std::make_shared< Variant >( refSeq, Region( "20", 56258612, 56258612 ), "GTGTGTGT" );
    const auto la = ins->getLeftAligned( 56258600 );
    const auto ra = ins->getRightAligned( 56258650 );

    BOOST_CHECK_EQUAL( la->region(), Region( "20", 56258610, 56258610 ) );
    BOOST_CHECK_EQUAL( ra->region(), Region( "20", 56258640, 56258640 ) );
}

BOOST_AUTO_TEST_CASE( insertionnWillRightAlignSingleRepeatUnit )
{
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 100, 110 ), "AAAAATTTTT" );
    wecall::variant::varPtr_t del = std::make_shared< Variant >( refSeq, Region( "1", 100, 100 ), "A" );

    auto del2 = del->getRightAligned( 109 );

    BOOST_CHECK_EQUAL( del2->contig(), "1" );
    BOOST_CHECK_EQUAL( del2->start(), 105 );
    BOOST_CHECK_EQUAL( del2->end(), 105 );
    BOOST_CHECK_EQUAL( del2->sequence(), wecall::utils::BasePairSequence( "A" ) );
    BOOST_CHECK( not del2->isFullyLeftAligned() );
}

BOOST_AUTO_TEST_CASE( insertionWillRightAlignMultiRepeatUnit )
{
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 100, 110 ), "ATATATACCC" );
    wecall::variant::varPtr_t var = std::make_shared< Variant >( refSeq, Region( "1", 100, 100 ), "AT" );

    auto var2 = var->getRightAligned( 109 );

    BOOST_CHECK_EQUAL( var2->contig(), "1" );
    BOOST_CHECK_EQUAL( var2->start(), 107 );
    BOOST_CHECK_EQUAL( var2->end(), 107 );
    BOOST_CHECK_EQUAL( var2->sequence(), wecall::utils::BasePairSequence( "TA" ) );
    BOOST_CHECK( not var2->isFullyLeftAligned() );

    BOOST_CHECK_EQUAL( var2->getStartEndRegions( var2->start() - 1, var2->end() + 1 ),
                       SetRegions( Region( "1", 106, 107 ) ) );
    BOOST_CHECK_EQUAL( var2->getStartEndRegions( 100, 110 ), SetRegions( Region( "1", 100, 107 ) ) );
}

BOOST_AUTO_TEST_CASE( testInsertionDoesJoinWithSNP )
{
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 10, 20 ), "CAGATTACAG" );

    const auto ins = std::make_shared< Variant >( refSeq, Region( "1", 11, 11 ), "A" );
    const auto snp = std::make_shared< Variant >( refSeq, Region( "1", 11, 12 ), "T" );

    const auto join = ins->join( snp );
    BOOST_REQUIRE( ins->joinable( snp ) );
    BOOST_CHECK_EQUAL( join->refSequence(), refSeq->subseq( Region( "1", 12, 12 ) ) );
    BOOST_CHECK_EQUAL( join->sequence(), "T" );
}

BOOST_AUTO_TEST_CASE( testSNPDoesJoinWithInsertion )
{
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 10, 20 ), "CAGATTACAG" );

    const auto snp = std::make_shared< Variant >( refSeq, Region( "1", 11, 12 ), "T" );
    const auto ins = std::make_shared< Variant >( refSeq, Region( "1", 12, 12 ), "A" );

    const auto join = snp->join( ins );
    BOOST_REQUIRE( snp->joinable( ins ) );
    BOOST_CHECK_EQUAL( join->refSequence(), refSeq->subseq( Region( "1", 11, 11 ) ) );
    BOOST_CHECK_EQUAL( join->sequence(), "T" );
}

BOOST_AUTO_TEST_CASE( testInsertionDoesJoinWithMNP )
{
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 10, 20 ), "CAGATTACAG" );

    const auto ins1 = std::make_shared< Variant >( refSeq, Region( "1", 11, 11 ), "A" );
    const auto snp = std::make_shared< Variant >( refSeq, Region( "1", 11, 13 ), "TT" );

    BOOST_CHECK( ins1->joinable( snp ) );
}

BOOST_AUTO_TEST_CASE( testInsertionDoesJoinWithDeletion )
{
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 10, 20 ), "CAGATTACAG" );

    const auto ins = std::make_shared< Variant >( refSeq, Region( "1", 11, 11 ), "C" );
    const auto del = std::make_shared< Variant >( refSeq, Region( "1", 11, 12 ), "" );

    auto join = ins->join( del );
    BOOST_REQUIRE( ins->joinable( del ) );
    BOOST_CHECK_EQUAL( join->refSequence(), refSeq->subseq( Region( "1", 11, 12 ) ) );
    BOOST_CHECK_EQUAL( join->sequence(), "C" );
}
