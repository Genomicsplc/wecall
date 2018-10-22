// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "variant/type/variant.hpp"
#include "io/read.hpp"
#include "alignment/cigar.hpp"
#include "caller/region.hpp"

using namespace wecall::variant;
using namespace wecall::alignment;
using namespace wecall::caller;
using wecall::utils::ReferenceSequence;

BOOST_AUTO_TEST_CASE( should_not_alter_untrimmable_variant )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 12 ), "AT" );
    auto var = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "CG" );
    auto trimmed = var->getTrimmed();

    BOOST_CHECK_EQUAL( *var, *trimmed );
}

BOOST_AUTO_TEST_CASE( should_trim_right_matching_sequence )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 12 ), "AT" );
    auto var = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "CT" );
    auto trimmed = var->getTrimmed();

    BOOST_CHECK_EQUAL( trimmed->refSequence(), referenceSequence->subseq( Region( "1", 10, 11 ) ) );
    BOOST_CHECK_EQUAL( trimmed->sequence(), "C" );
}

BOOST_AUTO_TEST_CASE( should_trim_left_matching_sequence )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 12 ), "TA" );
    auto var = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "TC" );
    auto trimmed = var->getTrimmed();

    BOOST_CHECK_EQUAL( trimmed->refSequence(), referenceSequence->subseq( Region( "1", 11, 12 ) ) );
    BOOST_CHECK_EQUAL( trimmed->sequence(), "C" );
}

BOOST_AUTO_TEST_CASE( should_trim_right_matching_sequence_before_left )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 14 ), "ATAT" );
    auto var = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "AT" );
    auto trimmed = var->getTrimmed();

    BOOST_CHECK_EQUAL( trimmed->refSequence(), referenceSequence->subseq( Region( "1", 10, 12 ) ) );
    BOOST_CHECK_EQUAL( trimmed->sequence(), "" );
}

BOOST_AUTO_TEST_CASE( should_trim_right_matching_sequence_then_remainder_of_left )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 15 ), "ATCAT" );
    auto var = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "AAT" );
    auto trimmed = var->getTrimmed();

    BOOST_CHECK_EQUAL( trimmed->refSequence(), referenceSequence->subseq( Region( "1", 11, 13 ) ) );
    BOOST_CHECK_EQUAL( trimmed->sequence(), "" );
}

BOOST_AUTO_TEST_CASE( should_get_correct_start_and_end_for_snp )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 11 ), "A" );
    auto snp = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "T" );
    BOOST_CHECK_EQUAL( snp->start(), 10 );
    BOOST_CHECK_EQUAL( snp->end(), 11 );
}

BOOST_AUTO_TEST_CASE( should_get_correct_start_and_end_for_mnp )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 2, 5 ), "BAA" );
    auto mnp = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "TCG" );
    BOOST_CHECK_EQUAL( mnp->start(), 2 );
    BOOST_CHECK_EQUAL( mnp->end(), 5 );
}

BOOST_AUTO_TEST_CASE( should_get_correct_start_and_end_for_ins )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 10 ), "" );
    auto ins = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "ATTT" );
    BOOST_CHECK_EQUAL( ins->start(), 10 );
    BOOST_CHECK_EQUAL( ins->end(), 10 );
}

BOOST_AUTO_TEST_CASE( should_get_correct_start_and_end_for_del )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 12 ), "AA" );
    auto del = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "" );
    BOOST_CHECK_EQUAL( del->start(), 10 );
    BOOST_CHECK_EQUAL( del->end(), 12 );
}

BOOST_AUTO_TEST_CASE( test_should_overlap_snp_with_itself )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 11 ), "A" );
    auto snp = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "T" );
    BOOST_CHECK( snp->overlaps( snp ) );
}

BOOST_AUTO_TEST_CASE( test_should_overlap_two_snps_at_same_position )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 11 ), "A" );
    auto snp1 = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "T" );
    auto snp2 = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "C" );

    BOOST_CHECK( snp1->overlaps( snp2 ) );
    BOOST_CHECK( snp2->overlaps( snp1 ) );
}

BOOST_AUTO_TEST_CASE( test_should_not_overlap_snps_at_different_positions )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 12 ), "AA" );
    auto snp1 = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 11 ), "T" );
    auto snp2 = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "T" );

    BOOST_CHECK( not snp1->overlaps( snp2 ) );
    BOOST_CHECK( not snp2->overlaps( snp1 ) );
}

BOOST_AUTO_TEST_CASE( test_should_not_overlap_snps_at_different_contigs )
{
    const auto referenceSequence1 = std::make_shared< ReferenceSequence >( Region( "1", 10, 11 ), "A" );
    const auto referenceSequence2 = std::make_shared< ReferenceSequence >( Region( "2", 10, 11 ), "A" );
    auto snp1 = std::make_shared< Variant >( referenceSequence1, referenceSequence1->region(), "T" );
    auto snp2 = std::make_shared< Variant >( referenceSequence2, referenceSequence2->region(), "G" );

    BOOST_CHECK( not snp1->overlaps( snp2 ) );
    BOOST_CHECK( not snp2->overlaps( snp1 ) );
}

BOOST_AUTO_TEST_CASE( test_should_overlap_two_mnp_and_snp_at_overlapping_position )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 12 ), "AA" );
    auto mnp = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 12 ), "TC" );
    auto snp = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "G" );

    BOOST_CHECK( mnp->overlaps( snp ) );
    BOOST_CHECK( snp->overlaps( mnp ) );
}

BOOST_AUTO_TEST_CASE( test_should_not_overlap_mnp_and_snp_at_different_positions )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 13 ), "AAT" );
    auto mnp = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 12 ), "TC" );
    auto snp = std::make_shared< Variant >( referenceSequence, Region( "1", 12, 13 ), "G" );

    BOOST_CHECK( not mnp->overlaps( snp ) );
    BOOST_CHECK( not snp->overlaps( mnp ) );
}

BOOST_AUTO_TEST_CASE( test_should_overlap_deletion_and_snp_at_same_positions )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 11 ), "A" );
    auto del = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "" );
    auto snp = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "G" );

    BOOST_CHECK( del->overlaps( snp ) );
    BOOST_CHECK( snp->overlaps( del ) );
}

BOOST_AUTO_TEST_CASE( test_should_not_overlap_deletion_and_snp_at_different_positions )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 13 ), "AAT" );
    auto del = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 12 ), "" );
    auto snp = std::make_shared< Variant >( referenceSequence, Region( "1", 12, 13 ), "G" );

    BOOST_CHECK( not del->overlaps( snp ) );
    BOOST_CHECK( not snp->overlaps( del ) );
}

BOOST_AUTO_TEST_CASE( test_should_not_overlap_insertion_and_snp )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 11 ), "A" );
    auto ins = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 10 ), "TT" );
    auto snp = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 11 ), "G" );

    BOOST_CHECK( not ins->overlaps( snp ) );
    BOOST_CHECK( not snp->overlaps( ins ) );
}

BOOST_AUTO_TEST_CASE( test_should_not_overlap_insertion_and_mnp )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 12 ), "AA" );
    auto ins = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 10 ), "TT" );
    auto mnp = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 12 ), "GG" );

    BOOST_CHECK( not ins->overlaps( mnp ) );
    BOOST_CHECK( not mnp->overlaps( ins ) );
}

BOOST_AUTO_TEST_CASE( test_should_not_overlap_insertion_and_deletion_at_same_position )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 11, 13 ), "AA" );
    auto ins = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 11 ), "TT" );
    auto del = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 13 ), "" );

    BOOST_CHECK( not ins->overlaps( del ) );
    BOOST_CHECK( not del->overlaps( ins ) );
}

BOOST_AUTO_TEST_CASE( test_should_overlap_insertion_and_deletion_at_overlapping_position )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 12 ), "AA" );
    auto ins = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 11 ), "TT" );
    auto del = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 12 ), "" );

    BOOST_CHECK( ins->overlaps( del ) );
    BOOST_CHECK( del->overlaps( ins ) );
}

BOOST_AUTO_TEST_CASE( should_not_overlap_insertion_before_region )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 15, 15 ), "" );
    auto ins = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "AAAAA" );
    BOOST_CHECK( not ins->overlaps( Region( "1", 15, 20 ) ) );
}

BOOST_AUTO_TEST_CASE( should_overlap_insertion_at_region_start )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 16, 16 ), "" );
    auto ins = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "AAAAA" );
    BOOST_CHECK( ins->overlaps( Region( "1", 15, 20 ) ) );
}

BOOST_AUTO_TEST_CASE( should_not_overlap_insertion_after_region )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 16, 16 ), "" );
    auto ins = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "AAAAA" );
    BOOST_CHECK( not ins->overlaps( Region( "1", 16, 20 ) ) );
}

BOOST_AUTO_TEST_CASE( should_overlap_deletion_touching_region_start )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 11, 16 ), "AAAAA" );
    auto del = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "" );
    BOOST_CHECK( del->overlaps( Region( "1", 15, 20 ) ) );
}

BOOST_AUTO_TEST_CASE( should_not_overlap_deletion_one_before_region_start )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 11, 16 ), "AAAAA" );
    auto del = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "" );
    BOOST_CHECK( not del->overlaps( Region( "1", 16, 20 ) ) );
}

BOOST_AUTO_TEST_CASE( should_overlap_deletion_touching_region_end )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 11, 16 ), "AAAAA" );
    auto del = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "" );
    BOOST_CHECK( del->overlaps( Region( "1", 0, 12 ) ) );
}

BOOST_AUTO_TEST_CASE( should_not_overlap_deletion_one_after_region_end )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 11, 16 ), "AAAAA" );
    auto del = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "" );
    BOOST_CHECK( not del->overlaps( Region( "1", 0, 11 ) ) );
}

BOOST_AUTO_TEST_CASE( should_not_overlap_mnp_with_region_on_left )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 11, 16 ), "AAAAA" );
    auto mnp = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "TAAAT" );
    BOOST_CHECK( not mnp->overlaps( Region( "1", 1, 11 ) ) );
}

BOOST_AUTO_TEST_CASE( should_overlap_mnp_with_region_on_left )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 11, 16 ), "AAAAA" );
    auto mnp = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "TAAAT" );
    BOOST_CHECK( mnp->overlaps( Region( "1", 1, 12 ) ) );
}

BOOST_AUTO_TEST_CASE( should_not_overlap_mnp_with_region_on_right )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 11, 16 ), "AAAAA" );
    auto mnp = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "TAAAT" );
    BOOST_CHECK( not mnp->overlaps( Region( "1", 16, 20 ) ) );
}

BOOST_AUTO_TEST_CASE( should_overlap_mnp_with_region_on_right )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 11, 16 ), "AAAAA" );
    auto mnp = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "TAAAT" );
    BOOST_CHECK( mnp->overlaps( Region( "1", 15, 20 ) ) );
}

BOOST_AUTO_TEST_CASE( should_not_overlap_del_before_with_read_with_matches )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 20 ), "AAAAAAAAAA" );
    // Variant deletes pos < 15.
    auto del = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 15 ), "" );
    wecall::io::Read read( std::string( 5, 'A' ), std::string( 5, 'Q' ), "", Cigar( "5M" ), 0, 15, 0, 0, 0, 0, 0,
                            referenceSequence );
    BOOST_CHECK( not del->interval().overlaps( read.getMaximalReadInterval() ) );
}

BOOST_AUTO_TEST_CASE( should_overlap_del_touching_start_of_read_with_matches )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 11, 21 ), "AAAAAAAAAA" );
    // Variant deletes pos 15.
    auto del = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 16 ), "" );
    wecall::io::Read read( std::string( 5, 'A' ), std::string( 5, 'Q' ), "", Cigar( "5M" ), 0, 15, 0, 0, 0, 0, 0,
                            referenceSequence );
    BOOST_CHECK( del->interval().overlaps( read.getMaximalReadInterval() ) );
}

BOOST_AUTO_TEST_CASE( should_not_overlap_del_before_with_read_with_insertion_at_start )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 15 ), "AAAAA" );
    // Variant deletes pos < 15.
    auto del = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "" );
    wecall::io::Read read( std::string( 5, 'A' ), std::string( 5, 'Q' ), "", Cigar( "1I4M" ), 0, 16, 0, 0, 0, 0, 0,
                            referenceSequence );
    BOOST_CHECK( not del->interval().overlaps( read.getMaximalReadInterval() ) );
}

BOOST_AUTO_TEST_CASE( should_overlap_del_touching_start_of_read_with_insertion_at_start )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 11, 16 ), "AAAAA" );
    // Variant deletes pos 15.
    auto del = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "" );
    // Read starts at aligned pos 16 BUT has an insertion before first piece of matching sequence in effect could be
    // flattened to start at 15.
    wecall::io::Read read( std::string( 5, 'A' ), std::string( 5, 'Q' ), "", Cigar( "1I4M" ), 0, 16, 0, 0, 0, 0, 0,
                            referenceSequence );
    BOOST_CHECK( del->interval().overlaps( read.getMaximalReadInterval() ) );
}

BOOST_AUTO_TEST_CASE( should_not_overlap_del_after_with_read_with_matches )
{
    const auto referenceSequence =
        std::make_shared< ReferenceSequence >( Region( "1", 14, 25 ), std::string( 25 - 14, 'A' ) );
    // Variant deletes pos >= 20.
    // Read last pos it 19.
    auto del = std::make_shared< Variant >( referenceSequence, Region( "1", 20, 25 ), "" );
    wecall::io::Read read( std::string( 5, 'A' ), std::string( 5, 'Q' ), "", Cigar( "5M" ), 0, 15, 0, 0, 0, 0, 0,
                            referenceSequence );
    BOOST_CHECK( not del->interval().overlaps( read.getMaximalReadInterval() ) );
}

BOOST_AUTO_TEST_CASE( should_overlap_del_touching_end_of_read_with_matches )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 14, 24 ), "AAAAAAAAAA" );
    // Variant deletes pos 19.
    // Read last pos is 19
    auto del = std::make_shared< Variant >( referenceSequence, Region( "1", 19, 24 ), "" );
    wecall::io::Read read( std::string( 5, 'A' ), std::string( 5, 'Q' ), "", Cigar( "5M" ), 0, 15, 0, 0, 0, 0, 0,
                            referenceSequence );
    BOOST_CHECK( del->interval().overlaps( read.getMaximalReadInterval() ) );
}

BOOST_AUTO_TEST_CASE( should_not_overlap_del_after_with_read_with_insertion_at_end )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 20, 25 ), "AAAAA" );
    // Variant deletes pos >= 20.
    // Read last aligned pos is 18 with one insertion.
    auto del = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "" );
    wecall::io::Read read( std::string( 5, 'A' ), std::string( 5, 'Q' ), "", Cigar( "4M1I" ), 0, 15, 0, 0, 0, 0, 0,
                            referenceSequence );
    BOOST_CHECK( not del->interval().overlaps( read.getMaximalReadInterval() ) );
}

BOOST_AUTO_TEST_CASE( should_overlap_del_touching_end_of_read_with_insertion_at_end )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 19, 24 ), "AAAAA" );
    // Variant deletes pos 19.
    // Read last aligned pos is 18 but has one bit of insertion meaning they "touch"
    auto del = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "" );
    wecall::io::Read read( std::string( 5, 'A' ), std::string( 5, 'Q' ), "", Cigar( "4M1I" ), 0, 15, 0, 0, 0, 0, 0,
                            referenceSequence );
    BOOST_CHECK( del->interval().overlaps( read.getMaximalReadInterval() ) );
}

BOOST_AUTO_TEST_CASE( shouldSortVariantsFirstOnContig )
{
    const auto referenceA = std::make_shared< ReferenceSequence >( Region( "A", 10, 11 ), "A" );
    const auto referenceB = std::make_shared< ReferenceSequence >( Region( "B", 10, 11 ), "A" );
    auto snpB = std::make_shared< Variant >( referenceB, referenceB->region(), "T" );
    auto snpA = std::make_shared< Variant >( referenceA, referenceA->region(), "T" );

    std::vector< varPtr_t > variants = {snpB, snpA};

    std::sort( variants.begin(), variants.end(), wecall::variant::varPtrComp() );

    std::vector< varPtr_t > expectedResult = {snpA, snpB};

    // Then
    BOOST_CHECK_EQUAL_COLLECTIONS( variants.begin(), variants.end(), expectedResult.begin(), expectedResult.end() );
}

BOOST_AUTO_TEST_CASE( shouldSortVariantsSecondOnLastRefPosBeforeSequence )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 13 ), "AAA" );
    auto del = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 13 ), "" );
    auto mnp = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 13 ), "TTT" );

    BOOST_CHECK_EQUAL( del->contig(), mnp->contig() );

    std::vector< varPtr_t > variants = {del, mnp};

    std::sort( variants.begin(), variants.end(), wecall::variant::varPtrComp() );

    std::vector< varPtr_t > expectedResult = {mnp, del};

    // Then
    BOOST_CHECK_EQUAL_COLLECTIONS( variants.begin(), variants.end(), expectedResult.begin(), expectedResult.end() );
}

BOOST_AUTO_TEST_CASE( shouldSortThirdlyEndPosition )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 13 ), "AAA" );
    auto del = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 12 ), "" );
    auto mnp = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 13 ), "TTT" );

    BOOST_CHECK_EQUAL( del->contig(), mnp->contig() );
    BOOST_CHECK_EQUAL( del->start(), mnp->start() );

    std::vector< varPtr_t > variants = {mnp, del};

    std::sort( variants.begin(), variants.end(), wecall::variant::varPtrComp() );

    std::vector< varPtr_t > expectedResult = {del, mnp};

    // Then
    BOOST_CHECK_EQUAL_COLLECTIONS( variants.begin(), variants.end(), expectedResult.begin(), expectedResult.end() );
}

BOOST_AUTO_TEST_CASE( shouldSortForthlyOnend )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 13 ), "AAA" );
    auto snp = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 11 ), "T" );
    auto mnp = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 13 ), "TTT" );

    BOOST_CHECK_EQUAL( snp->contig(), mnp->contig() );
    BOOST_CHECK_EQUAL( snp->start(), mnp->start() );
    BOOST_CHECK_EQUAL( snp->zeroIndexedVcfPosition(), mnp->zeroIndexedVcfPosition() );

    std::vector< varPtr_t > variants = {mnp, snp};

    std::sort( variants.begin(), variants.end(), wecall::variant::varPtrComp() );

    std::vector< varPtr_t > expectedResult = {snp, mnp};

    // Then
    BOOST_CHECK_EQUAL_COLLECTIONS( variants.begin(), variants.end(), expectedResult.begin(), expectedResult.end() );
}

BOOST_AUTO_TEST_CASE( shouldSortFiftlyOnSequenceLength )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 11, 11 ), "" );
    auto ins_1 = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "A" );
    auto ins_2 = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "AAA" );

    BOOST_CHECK_EQUAL( ins_1->contig(), ins_2->contig() );
    BOOST_CHECK_EQUAL( ins_1->start(), ins_2->start() );
    BOOST_CHECK_EQUAL( ins_1->zeroIndexedVcfPosition(), ins_2->zeroIndexedVcfPosition() );
    BOOST_CHECK_EQUAL( ins_1->end(), ins_2->end() );

    std::vector< varPtr_t > variants = {ins_2, ins_1};

    std::sort( variants.begin(), variants.end(), wecall::variant::varPtrComp() );

    std::vector< varPtr_t > expectedResult = {ins_1, ins_2};

    // Then
    BOOST_CHECK_EQUAL_COLLECTIONS( variants.begin(), variants.end(), expectedResult.begin(), expectedResult.end() );
}

BOOST_AUTO_TEST_CASE( shouldSortFinallyOnSequence )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 11, 11 ), "" );
    auto mnp_1 = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "TAT" );
    auto mnp_2 = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "TTT" );

    BOOST_CHECK_EQUAL( mnp_1->contig(), mnp_2->contig() );
    BOOST_CHECK_EQUAL( mnp_1->start(), mnp_2->start() );
    BOOST_CHECK_EQUAL( mnp_1->start(), mnp_2->start() );
    BOOST_CHECK_EQUAL( mnp_1->zeroIndexedVcfPosition(), mnp_2->zeroIndexedVcfPosition() );
    BOOST_CHECK_EQUAL( mnp_1->end(), mnp_2->end() );
    BOOST_CHECK_EQUAL( mnp_1->sequenceLength(), mnp_2->sequenceLength() );

    std::vector< varPtr_t > variants = {mnp_2, mnp_1};

    std::sort( variants.begin(), variants.end(), wecall::variant::varPtrComp() );

    std::vector< varPtr_t > expectedResult = {mnp_1, mnp_2};

    // Then
    BOOST_CHECK_EQUAL_COLLECTIONS( variants.begin(), variants.end(), expectedResult.begin(), expectedResult.end() );
}

BOOST_AUTO_TEST_CASE( shouldSortCollectionAtSameVCFPosition )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 10, 13 ), "AAA" );
    auto mnp = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 13 ), "TAT" );
    auto snp = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 11 ), "C" );
    auto del = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "" );
    auto ins = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 11 ), "A" );

    BOOST_CHECK_EQUAL( mnp->zeroIndexedVcfPosition(), 10 );
    BOOST_CHECK_EQUAL( snp->zeroIndexedVcfPosition(), 10 );
    BOOST_CHECK_EQUAL( del->zeroIndexedVcfPosition(), 10 );
    BOOST_CHECK_EQUAL( ins->zeroIndexedVcfPosition(), 10 );

    std::vector< varPtr_t > variants = {del, ins, mnp, snp};

    std::sort( variants.begin(), variants.end(), wecall::variant::varPtrComp() );

    std::vector< varPtr_t > expectedResult = {snp, mnp, ins, del};

    // Then
    BOOST_CHECK_EQUAL_COLLECTIONS( variants.begin(), variants.end(), expectedResult.begin(), expectedResult.end() );
}

BOOST_AUTO_TEST_CASE( shouldSortCollectionWithLastRefBeforeSequence )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 11, 14 ), "AAA" );
    auto mnp = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 14 ), "TAT" );
    auto snp = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "C" );
    auto del = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 12 ), "" );
    auto ins = std::make_shared< Variant >( referenceSequence, Region( "1", 11, 11 ), "A" );

    BOOST_CHECK_EQUAL( mnp->start(), 11 );
    BOOST_CHECK_EQUAL( snp->start(), 11 );
    BOOST_CHECK_EQUAL( del->start(), 11 );
    BOOST_CHECK_EQUAL( ins->start(), 11 );

    std::vector< varPtr_t > variants = {del, ins, mnp, snp};

    std::sort( variants.begin(), variants.end(), wecall::variant::varPtrComp() );

    std::vector< varPtr_t > expectedResult = {ins, del, snp, mnp};

    // Then
    BOOST_CHECK_EQUAL_COLLECTIONS( variants.begin(), variants.end(), expectedResult.begin(), expectedResult.end() );
}

BOOST_AUTO_TEST_CASE( testSNPVariantSplit )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "contig", 1, 2 ), "A" );
    Variant snp( referenceSequence, referenceSequence->region(), "T" );

    const auto splat = snp.split();
    BOOST_REQUIRE_EQUAL( splat.size(), 1 );

    BOOST_CHECK_EQUAL( splat[0]->refSequence(), referenceSequence->subseq( Region( "contig", 1, 2 ) ) );
    BOOST_CHECK_EQUAL( splat[0]->sequence(), "T" );
}

BOOST_AUTO_TEST_CASE( testPureDeletionVariantSplit )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 8 ), "AAAAAAAA" );
    Variant var( referenceSequence, referenceSequence->region(), "" );

    auto splat = var.split();
    BOOST_REQUIRE_EQUAL( splat.size(), 1 );

    BOOST_CHECK_EQUAL( *splat[0], var );
}

BOOST_AUTO_TEST_CASE( testPureInsertionVariantSplit )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 0 ), "" );
    Variant var( referenceSequence, referenceSequence->region(), "TA" );

    auto splat = var.split();
    BOOST_REQUIRE_EQUAL( splat.size(), 1 );

    BOOST_CHECK_EQUAL( *splat[0], var );
}

BOOST_AUTO_TEST_CASE( testMNPVariantSplit )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "contig", 1, 3 ), "AT" );
    Variant mnp( referenceSequence, referenceSequence->region(), "TC" );

    const auto splat = mnp.split();

    BOOST_REQUIRE_EQUAL( splat.size(), 2 );

    BOOST_CHECK_EQUAL( *splat[0], Variant( referenceSequence, Region( "contig", 1, 2 ), "T" ) );
    BOOST_CHECK_EQUAL( *splat[1], Variant( referenceSequence, Region( "contig", 2, 3 ), "C" ) );
}

BOOST_AUTO_TEST_CASE( testComplexDeletionVariantSplit )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    Variant complexDel( referenceSequence, referenceSequence->region(), "TC" );

    auto splat = complexDel.split();
    BOOST_REQUIRE_EQUAL( splat.size(), 3 );

    BOOST_CHECK_EQUAL( splat[0]->refSequence(), referenceSequence->subseq( Region( "1", 0, 8 ) ) );
    BOOST_CHECK_EQUAL( splat[0]->sequence(), "" );

    BOOST_CHECK_EQUAL( splat[1]->refSequence(), referenceSequence->subseq( Region( "1", 8, 9 ) ) );
    BOOST_CHECK_EQUAL( splat[1]->sequence(), "T" );

    BOOST_CHECK_EQUAL( splat[2]->refSequence(), referenceSequence->subseq( Region( "1", 9, 10 ) ) );
    BOOST_CHECK_EQUAL( splat[2]->sequence(), "C" );
}

BOOST_AUTO_TEST_CASE( testComplexInsertionVariantSplit )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 2 ), "AA" );
    Variant complexIns( referenceSequence, referenceSequence->region(), "TAAC" );

    auto splat = complexIns.split();
    BOOST_REQUIRE_EQUAL( splat.size(), 2 );
    BOOST_CHECK_EQUAL( splat[0]->refSequence(), referenceSequence->subseq( Region( "1", 0, 0 ) ) );
    BOOST_CHECK_EQUAL( splat[0]->sequence(), "TA" );

    BOOST_CHECK_EQUAL( splat[1]->refSequence(), referenceSequence->subseq( Region( "1", 1, 2 ) ) );
    BOOST_CHECK_EQUAL( splat[1]->sequence(), "C" );
}

BOOST_AUTO_TEST_CASE( testComplexDeletionVariantSplit2 )
{
    const auto referenceSequence =
        std::make_shared< ReferenceSequence >( Region( "19", 41611054, 41611067 ), "AAAACACCAAAAG" );
    Variant complexDel( referenceSequence, referenceSequence->region(), "CTTCCCTCT" );

    auto splat = complexDel.split();

    BOOST_REQUIRE_EQUAL( splat.size(), 8 );

    BOOST_CHECK_EQUAL( splat[0]->refSequence(), referenceSequence->subseq( Region( "19", 41611054, 41611058 ) ) );
    BOOST_CHECK_EQUAL( splat[0]->sequence(), "" );

    BOOST_CHECK_EQUAL( splat[1]->refSequence(), referenceSequence->subseq( Region( "19", 41611059, 41611060 ) ) );
    BOOST_CHECK_EQUAL( splat[1]->sequence(), "T" );

    BOOST_CHECK_EQUAL( splat[2]->refSequence(), referenceSequence->subseq( Region( "19", 41611060, 41611061 ) ) );
    BOOST_CHECK_EQUAL( splat[2]->sequence(), "T" );

    BOOST_CHECK_EQUAL( splat[3]->refSequence(), referenceSequence->subseq( Region( "19", 41611062, 41611063 ) ) );
    BOOST_CHECK_EQUAL( splat[3]->sequence(), "C" );

    BOOST_CHECK_EQUAL( splat[4]->refSequence(), referenceSequence->subseq( Region( "19", 41611063, 41611064 ) ) );
    BOOST_CHECK_EQUAL( splat[4]->sequence(), "C" );

    BOOST_CHECK_EQUAL( splat[5]->refSequence(), referenceSequence->subseq( Region( "19", 41611064, 41611065 ) ) );
    BOOST_CHECK_EQUAL( splat[5]->sequence(), "T" );

    BOOST_CHECK_EQUAL( splat[6]->refSequence(), referenceSequence->subseq( Region( "19", 41611065, 41611066 ) ) );
    BOOST_CHECK_EQUAL( splat[6]->sequence(), "C" );

    BOOST_CHECK_EQUAL( splat[7]->refSequence(), referenceSequence->subseq( Region( "19", 41611066, 41611067 ) ) );
    BOOST_CHECK_EQUAL( splat[7]->sequence(), "T" );
}

BOOST_AUTO_TEST_CASE( testFromBreakpointToggle )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    Variant variant( referenceSequence, referenceSequence->region(), "TC" );
    BOOST_CHECK( not variant.isFromBreakpoint() );
    variant.setFromBreakpoint();
    BOOST_CHECK( variant.isFromBreakpoint() );
}

BOOST_AUTO_TEST_CASE( testVariantOnlyRemovableIfContainedInOthersRegion )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    const auto variant1 = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "TC" );
    const auto variant2 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 10 ), "TC" );

    BOOST_CHECK( variant1->removable( variant1 ) );
    BOOST_CHECK( variant1->removable( variant2 ) );

    const auto expectEmpty = variant1->remove( variant1 );
    const auto expectPureDeletion = variant1->remove( variant2 );

    BOOST_CHECK( expectEmpty->empty() );
    BOOST_CHECK_EQUAL( *expectPureDeletion, Variant( referenceSequence, Region( "1", 0, 1 ), "" ) );

    BOOST_CHECK( not variant2->removable( variant1 ) );
    BOOST_CHECK_THROW( variant2->remove( variant1 ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( testVariantOnlyRemovableIfStartsOrEndsAtSameLocation )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    const auto variant1 = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "TC" );
    const auto variant2 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 9 ), "TC" );

    BOOST_CHECK( not variant1->removable( variant2 ) );
    BOOST_CHECK_THROW( variant1->remove( variant2 ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( testVariantOnlyRemovableIfAlternateSequenceMatches )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    const auto variant1 = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "TC" );
    const auto variant2 = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 5 ), "T" );
    const auto variant3 = std::make_shared< Variant >( referenceSequence, Region( "1", 5, 10 ), "C" );
    const auto variant4 = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 5 ), "C" );
    const auto variant5 = std::make_shared< Variant >( referenceSequence, Region( "1", 5, 10 ), "T" );

    BOOST_CHECK( variant1->removable( variant2 ) );
    BOOST_CHECK( variant1->removable( variant3 ) );
    BOOST_CHECK( not variant1->removable( variant4 ) );
    BOOST_CHECK( not variant1->removable( variant5 ) );

    BOOST_CHECK_EQUAL( *variant1->remove( variant2 ), Variant( referenceSequence, Region( "1", 5, 10 ), "C" ) );
    BOOST_CHECK_EQUAL( *variant1->remove( variant3 ), Variant( referenceSequence, Region( "1", 0, 5 ), "T" ) );

    BOOST_CHECK_THROW( variant1->remove( variant4 ), wecall::utils::wecall_exception );
    BOOST_CHECK_THROW( variant1->remove( variant5 ), wecall::utils::wecall_exception );
}

BOOST_AUTO_TEST_CASE( testShouldTrimResultOfRemovingVariant )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    const auto variant1 = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "TAAAAAAAAC" );
    const auto variant2 = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "T" );
    const auto variant3 = std::make_shared< Variant >( referenceSequence, Region( "1", 9, 10 ), "C" );

    BOOST_CHECK_EQUAL( variant1->sequence().str(), "TAAAAAAAAC" );
    const auto withoutVariant2 = variant1->remove( variant2 );
    BOOST_CHECK_EQUAL( withoutVariant2->sequence().str(), "C" );
    const auto expectEmpty = withoutVariant2->remove( variant3 );
    BOOST_CHECK( expectEmpty->empty() );
}

BOOST_AUTO_TEST_CASE( testShouldCacheReadsInVariantClass )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 5, 17 ), "GGGGGATGGGGG" );
    auto var = std::make_shared< Variant >( referenceSequence, Region( "1", 10, 12 ), "CG" );
    auto readPtr1 = std::make_shared< wecall::io::Read >( "EDWARD", "ADRIAN", "", wecall::alignment::Cigar( "6M" ), 0,
                                                           10, 0, 0, 0, 0, 0, referenceSequence );
    auto readPtr2 = std::make_shared< wecall::io::Read >( "ADRIAN", "EDWARD", "", wecall::alignment::Cigar( "6M" ), 0,
                                                           10, 0, 0, 0, 0, 0, referenceSequence );

    var->addRead( readPtr1 );
    var->addRead( readPtr2 );

    const std::vector< wecall::io::readPtr_t > expected = {readPtr1, readPtr2};
    const auto actual = var->getReads();

    BOOST_CHECK_EQUAL_COLLECTIONS( expected.begin(), expected.end(), actual.begin(), actual.end() );
}
