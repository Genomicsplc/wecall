// All content Copyright (C) 2018 Genomics plc
#include "variant/type/variant.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <string>

using echidna::caller::Region;
using echidna::utils::ReferenceSequence;
using echidna::variant::Variant;
using echidna::caller::SetRegions;

BOOST_AUTO_TEST_CASE( testDeletion )
{
    using namespace echidna::variant;

    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 1000001, 1000005 ), "GGCA" );
    varPtr_t theDeletion( new Variant( referenceSequence, referenceSequence->region(), "", false ) );

    BOOST_CHECK_EQUAL( theDeletion->sequence(), std::string( "" ) );
    BOOST_CHECK_EQUAL( theDeletion->refSequence().sequence(), referenceSequence->sequence() );

    BOOST_CHECK_EQUAL( theDeletion->start(), 1000001 );
    BOOST_CHECK_EQUAL( theDeletion->end(), 1000005 );

    auto expectedPrior = 1e-4 * pow( 0.8, referenceSequence->size() );
    echidna::variant::setDefaultPriors( {theDeletion} );

    BOOST_CHECK_CLOSE( theDeletion->prior(), expectedPrior, 1e-5 );
}

BOOST_AUTO_TEST_CASE( leftAlignDel )
{
    // Test left-aligning deletions
    using namespace echidna::variant;

    // Reference genome 37, chrom 1 from 1000000 - 1000099
    //
    //                          0         1         2         3         4         5         6         7         8 9
    //                          0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    const auto refSeq = std::make_shared< ReferenceSequence >(
        Region( "1", 1000000, 1000100 ),
        "GGGCACAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACGGGAAGGAA" );

    varPtr_t del1( new Variant( refSeq, Region( "1", 1000053, 1000054 ), "", false ) );
    varPtr_t del2( new Variant( refSeq, Region( "1", 1000058, 1000059 ), "", false ) );
    varPtr_t del3( new Variant( refSeq, Region( "1", 1000058, 1000060 ), "", false ) );
    varPtr_t del4( new Variant( refSeq, Region( "1", 1000058, 1000061 ), "", false ) );
    varPtr_t del5( new Variant( refSeq, Region( "1", 1000058, 1000070 ), "", false ) );
    varPtr_t del6( new Variant( refSeq, Region( "1", 1000058, 1000068 ), "", false ) );
    varPtr_t del7( new Variant( refSeq, Region( "1", 1000055, 1000056 ), "", false ) );

    // Chec del1 left aligns correctly
    del1 = del1->getLeftAligned( 1000000 );
    BOOST_CHECK_EQUAL( del1->start(), 1000053 );
    BOOST_CHECK_EQUAL( del1->end(), 1000054 );
    BOOST_CHECK_EQUAL( del1->refSequence().sequence(), "C" );

    // Chec del2 left aligns correctly
    del2 = del2->getLeftAligned( 1000000 );
    BOOST_CHECK_EQUAL( del2->start(), 1000058 );
    BOOST_CHECK_EQUAL( del2->end(), 1000059 );
    BOOST_CHECK_EQUAL( del2->refSequence().sequence(), "T" );

    // Chec del3 left aligns correctly
    del3 = del3->getLeftAligned( 1000000 );
    BOOST_CHECK_EQUAL( del3->start(), 1000058 );
    BOOST_CHECK_EQUAL( del3->end(), 1000060 );
    BOOST_CHECK_EQUAL( del3->refSequence().sequence(), "TG" );

    // Chec del4 left aligns correctly
    del4 = del4->getLeftAligned( 1000000 );
    BOOST_CHECK_EQUAL( del4->start(), 1000057 );
    BOOST_CHECK_EQUAL( del4->end(), 1000060 );
    BOOST_CHECK_EQUAL( del4->refSequence().sequence(), "CTG" );

    // Chec del5 left aligns correctly
    del5 = del5->getLeftAligned( 1000000 );
    BOOST_CHECK_EQUAL( del5->start(), 1000057 );
    BOOST_CHECK_EQUAL( del5->end(), 1000069 );
    BOOST_CHECK_EQUAL( del5->refSequence().sequence(), "CTGCCAGGGAGA" );

    // Chec del6 left aligns correctly
    del6 = del6->getLeftAligned( 1000000 );
    BOOST_CHECK_EQUAL( del6->start(), 1000058 );
    BOOST_CHECK_EQUAL( del6->end(), 1000068 );
    BOOST_CHECK_EQUAL( del6->refSequence().sequence(), "TGCCAGGGAG" );

    // Chec del7 left aligns correctly
    del7 = del7->getLeftAligned( 1000000 );
    BOOST_CHECK_EQUAL( del7->start(), 1000053 );
    BOOST_CHECK_EQUAL( del7->end(), 1000054 );
    BOOST_CHECK_EQUAL( del7->refSequence().sequence(), "C" );
}

BOOST_AUTO_TEST_CASE( minPosForLeftAlignmentAtSamePosShouldNotAffectWhenNoAlignmentDoneDel )
{
    // Test left-aligning deletions
    using namespace echidna::variant;

    // Reference genome 37, chrom 1 from 1000000 - 1000099
    //
    //                          0         1         2         3         4         5         6         7         8 9
    //                          0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    const auto refSeq = std::make_shared< ReferenceSequence >(
        Region( "1", 1000000, 1000100 ),
        "GGGCACAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACGGGAAGGAA" );

    varPtr_t del = std::make_shared< Variant >( refSeq, Region( "1", 1000054, 1000055 ), "", false );

    // Chec ins1 left aligns correctly
    const auto del2 = del->getLeftAligned( 1000052 );
    BOOST_CHECK_EQUAL( del2->refSequence().sequence(), std::string( "C" ) );
    BOOST_CHECK_EQUAL( del2->start(), 1000053 );
    BOOST_CHECK_EQUAL( del2->end(), 1000054 );

    const auto del3 = del2->getRightAligned( 1000055 );

    BOOST_CHECK_EQUAL( del->sequenceLengthInRef(), del3->sequenceLengthInRef() );
    BOOST_CHECK_EQUAL( del->sequence(), del3->sequence() );

    echidna::variant::setDefaultPriors( {del, del3} );
    BOOST_CHECK_CLOSE( del->prior(), del3->prior(), 1e-5 );
}

BOOST_AUTO_TEST_CASE( deletionWontLeftAlignWhenBoundaryInsideRepetitiveRegion )
{
    // Test left-aligning deletions
    using namespace echidna::variant;

    // Reference genome 37, chrom 1 from 1000000 - 1000099
    //
    //                          0         1         2         3         4         5         6         7         8 9
    //                          0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    const auto refSeq = std::make_shared< ReferenceSequence >(
        Region( "1", 1000000, 1000100 ),
        "GGGCACAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACGGGAAGGAA" );

    varPtr_t del = std::make_shared< Variant >( refSeq, Region( "1", 1000057, 1000058 ), "", false );

    // Check del does not left align:  If block boundary is denoted by | , the situation is
    // AGA|CCCCCTG  which can't be left aligned -- can't move away from homopolymer.
    auto del2 = del->getLeftAligned( 1000053 );
    // Should still report the variant at the smallest possible position.
    BOOST_CHECK_EQUAL( del2->refSequence().sequence(), std::string( "C" ) );
    BOOST_CHECK_EQUAL( del2->start(), 1000053 );
    BOOST_CHECK_EQUAL( del2->end(), 1000054 );

    const auto del3 = del2->getRightAligned( del->end() );

    BOOST_CHECK_EQUAL( del->sequenceLengthInRef(), del3->sequenceLengthInRef() );
    BOOST_CHECK_EQUAL( del->sequence(), del3->sequence() );
    echidna::variant::setDefaultPriors( {del, del3} );
    BOOST_CHECK_CLOSE( del->prior(), del3->prior(), 1e-5 );
}

BOOST_AUTO_TEST_CASE( deletionWillLeftAlignWhenBoundaryInsideRepetitiveRegion )
{
    // Test left-aligning deletions
    using namespace echidna::variant;

    const auto refSeq = std::make_shared< ReferenceSequence >(
        Region( "1", 1000000, 1000100 ),
        "GGGCACAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACGGGAAGGAA" );

    varPtr_t del( new Variant( refSeq, Region( "1", 1000057, 1000058 ), "", false ) );

    // Check del does left align:  If block boundary is denoted by | , the situation is
    // AG|ACCCCCTG  which can be left aligned -- can move away from homopolymer.
    auto del2 = del->getLeftAligned( 1000052 );

    const auto del3 = del2->getRightAligned( del->end() );

    BOOST_CHECK_EQUAL( del->sequenceLengthInRef(), del3->sequenceLengthInRef() );
    BOOST_CHECK_EQUAL( del->sequence(), del3->sequence() );
    echidna::variant::setDefaultPriors( {del, del3} );
    BOOST_CHECK_CLOSE( del->prior(), del3->prior(), 1e-5 );
}

BOOST_AUTO_TEST_CASE( deletionWillRightAlignSingleRepeatUnit )
{
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 100, 112 ), "AAAAATTTTTAC" );
    echidna::variant::varPtr_t del = std::make_shared< Variant >( refSeq, Region( "1", 101, 102 ), "" );

    auto del2 = del->getRightAligned( 109 );

    BOOST_CHECK_EQUAL( del2->contig(), "1" );
    BOOST_CHECK_EQUAL( del2->start(), 104 );
    BOOST_CHECK_EQUAL( del2->end(), 105 );
    BOOST_CHECK_EQUAL( del2->refSequence().sequence(), echidna::utils::BasePairSequence( "A" ) );
    BOOST_CHECK( not del2->isFullyLeftAligned() );

    BOOST_CHECK_EQUAL( del2->getStartEndRegions( del2->start() - 1, del2->end() + 1 ),
                       SetRegions( Region( "1", 103, 105 ) ) );
    BOOST_CHECK_EQUAL( del2->getStartEndRegions( del2->start() - 4, del2->end() + 4 ),
                       SetRegions( Region( "1", 100, 105 ) ) );
    BOOST_CHECK_EQUAL( del2->getStartEndRegions( 100, 112 ), SetRegions( Region( "1", 100, 105 ) ) );
}

BOOST_AUTO_TEST_CASE( deletionWillRightAlignMultiRepeatUnit )
{
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 100, 110 ), "ATATATACCC" );
    echidna::variant::varPtr_t del = std::make_shared< Variant >( refSeq, Region( "1", 101, 103 ), "" );

    auto del2 = del->getRightAligned( 109 );

    BOOST_CHECK_EQUAL( del2->contig(), "1" );
    BOOST_CHECK_EQUAL( del2->start(), 105 );
    BOOST_CHECK_EQUAL( del2->end(), 107 );
    BOOST_CHECK_EQUAL( del2->refSequence().sequence(), echidna::utils::BasePairSequence( "TA" ) );
    BOOST_CHECK( not del2->isFullyLeftAligned() );
}

BOOST_AUTO_TEST_CASE( deletionWillLeftAlignInComplexRepeatRegion )
{
    const auto refSeq = std::make_shared< ReferenceSequence >(
        Region( "1", 51, 108 ), "TGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA" );

    using namespace echidna::variant;
    varPtr_t ins( new Variant( refSeq, Region( "1", 101, 108 ), "", false ) );
    //                                    0      1      2      3      4      5      6
    ins = ins->getLeftAligned( 51 );
    BOOST_CHECK_EQUAL( ins->refSequence().sequence(), std::string( "GATTACA" ) );
}

BOOST_AUTO_TEST_CASE( deletionWillLeftAlignInComplexRepeatRegionWithAllRefStringOfRepeatUnits )
{
    const auto refSeq = std::make_shared< ReferenceSequence >(
        Region( "1", 52, 108 ), "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA" );

    using namespace echidna::variant;
    varPtr_t del( new Variant( refSeq, Region( "1", 101, 108 ), "", false ) );
    //                                   0      1      2      3      4      5      6
    del = del->getLeftAligned( 52 );
    BOOST_CHECK_EQUAL( del->refSequence().sequence(), std::string( "GATTACA" ) );
    BOOST_CHECK( not del->isFullyLeftAligned() );
}

BOOST_AUTO_TEST_CASE( deletionWillLeftAlignInComplexRepeatRegionWithChangeToSequence )
{
    const auto refSeq = std::make_shared< ReferenceSequence >(
        Region( "1", 50, 108 ), "TAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA" );

    using namespace echidna::variant;
    varPtr_t ins( new Variant( refSeq, Region( "1", 101, 108 ), "", false ) );
    //                                     0      1      2      3      4      5      6
    const ReferenceSequence refStringBefore = refSeq->subseq( Region( "1", 50, 101 ) );

    ins = ins->getLeftAligned( 50 );
    BOOST_CHECK_EQUAL( ins->refSequence().sequence(), std::string( "AGATTAC" ) );
    BOOST_CHECK( ins->isFullyLeftAligned() );
}

BOOST_AUTO_TEST_CASE( deletionWillLeftAlignInComplexRepeatRegionWithChangeToSequenceAndAllRefStringMatchingRepeatUnits )
{
    const auto refSeq = std::make_shared< ReferenceSequence >(
        Region( "1", 50, 108 ), "CAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA" );
    using namespace echidna::variant;
    varPtr_t ins( new Variant( refSeq, Region( "1", 101, 108 ), "", false ) );
    //                                     0      1      2      3      4      5      6
    const ReferenceSequence refStringBefore = refSeq->subseq( Region( "1", 50, 101 ) );

    ins = ins->getLeftAligned( 50 );
    BOOST_CHECK_EQUAL( ins->refSequence().sequence(), std::string( "CAGATTA" ) );
}

BOOST_AUTO_TEST_CASE( testDeletionJoinsWithDeletionOnlyOnRight )
{
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 10, 20 ), "CAGATTACAG" );

    const auto del1 = std::make_shared< Variant >( refSeq, Region( "1", 10, 11 ), "" );
    const auto del2 = std::make_shared< Variant >( refSeq, Region( "1", 11, 12 ), "" );
    const auto del3 = std::make_shared< Variant >( refSeq, Region( "1", 13, 15 ), "" );

    const auto join = del1->join( del2 );

    BOOST_CHECK( del1->joinable( del2 ) );
    BOOST_CHECK_EQUAL( join->region(), Region( "1", 10, 12 ) );
    BOOST_CHECK_EQUAL( join->sequence(), "" );

    BOOST_CHECK( not del2->joinable( del1 ) );
    BOOST_CHECK( not del1->joinable( del3 ) );
}

BOOST_AUTO_TEST_CASE( testDeletionDoesJoinWithSNP )
{
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 10, 20 ), "CAGATTACAG" );

    const auto del = std::make_shared< Variant >( refSeq, Region( "1", 10, 11 ), "" );
    const auto snp = std::make_shared< Variant >( refSeq, Region( "1", 11, 12 ), "T" );

    const auto join = del->join( snp );

    BOOST_REQUIRE( del->joinable( snp ) );
    BOOST_CHECK_EQUAL( join->refSequence(), refSeq->subseq( Region( "1", 10, 12 ) ) );
    BOOST_CHECK_EQUAL( join->sequence(), "T" );
}

BOOST_AUTO_TEST_CASE( testSNPDoesJoinWithDeletion )
{
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 10, 20 ), "CAGATTACAG" );

    const auto snp = std::make_shared< Variant >( refSeq, Region( "1", 11, 12 ), "T" );
    const auto del = std::make_shared< Variant >( refSeq, Region( "1", 12, 13 ), "" );

    const auto join = snp->join( del );
    BOOST_REQUIRE( snp->joinable( del ) );
    BOOST_CHECK_EQUAL( join->refSequence(), refSeq->subseq( Region( "1", 11, 13 ) ) );
    BOOST_CHECK_EQUAL( join->sequence(), "T" );
}

BOOST_AUTO_TEST_CASE( testDeletionDoesJoinWithMNP )
{
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 10, 20 ), "CAGATTACAG" );

    const auto del1 = std::make_shared< Variant >( refSeq, Region( "1", 10, 11 ), "" );
    const auto mnp = std::make_shared< Variant >( refSeq, Region( "1", 11, 13 ), "TT" );

    BOOST_CHECK( del1->joinable( mnp ) );
}
