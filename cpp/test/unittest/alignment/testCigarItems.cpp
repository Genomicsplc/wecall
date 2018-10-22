// All content Copyright (C) 2018 Genomics plc
#include "unittest/vcf/VCFTestUtils.hpp"

#include "variant/type/variant.hpp"
#include "alignment/cigarItems.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <iostream>

using namespace wecall::alignment;
using namespace wecall::variant;
using wecall::caller::Region;
using wecall::utils::ReferenceSequence;

// Test get lengths from cigar items.

BOOST_AUTO_TEST_CASE( shouldGetLengthsForMatch )
{
    CigarMatch cigar( 2 );
    BOOST_CHECK_EQUAL( cigar.length(), 2 );
    BOOST_CHECK_EQUAL( cigar.lengthInRef(), 2 );
    BOOST_CHECK_EQUAL( cigar.lengthInSeq(), 2 );
}

BOOST_AUTO_TEST_CASE( shouldGetLengthsForDeletion )
{
    CigarDeletion cigar( 2 );
    BOOST_CHECK_EQUAL( cigar.length(), 2 );
    BOOST_CHECK_EQUAL( cigar.lengthInRef(), 2 );
    BOOST_CHECK_EQUAL( cigar.lengthInSeq(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldGetLengthsForInsertion )
{
    CigarInsertion cigar( 2 );
    BOOST_CHECK_EQUAL( cigar.length(), 2 );
    BOOST_CHECK_EQUAL( cigar.lengthInRef(), 0 );
    BOOST_CHECK_EQUAL( cigar.lengthInSeq(), 2 );
}

BOOST_AUTO_TEST_CASE( shouldGetLengthsForSkip )
{
    CigarSkip cigar( 2 );
    BOOST_CHECK_EQUAL( cigar.length(), 2 );
    BOOST_CHECK_EQUAL( cigar.lengthInRef(), 2 );
    BOOST_CHECK_EQUAL( cigar.lengthInSeq(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldGetLengthsForSoftClip )
{
    CigarSoftClip cigar( 2 );
    BOOST_CHECK_EQUAL( cigar.length(), 2 );
    BOOST_CHECK_EQUAL( cigar.lengthInRef(), 0 );
    BOOST_CHECK_EQUAL( cigar.lengthInSeq(), 2 );  // Keeps the sequence
}

BOOST_AUTO_TEST_CASE( shouldGetLengthsForHardClip )
{
    CigarHardClip cigar( 2 );
    BOOST_CHECK_EQUAL( cigar.length(), 2 );
    BOOST_CHECK_EQUAL( cigar.lengthInRef(), 0 );
    BOOST_CHECK_EQUAL( cigar.lengthInSeq(), 0 );  // Has also removed the sequence
}

BOOST_AUTO_TEST_CASE( shouldGetLengthsForPadding )
{
    CigarPadding cigar( 2 );
    BOOST_CHECK_EQUAL( cigar.length(), 2 );
    BOOST_CHECK_EQUAL( cigar.lengthInRef(), 0 );
    BOOST_CHECK_EQUAL( cigar.lengthInSeq(), 0 );
}

// Test get relative ref positions.

BOOST_AUTO_TEST_CASE( shouldGetTheCorrectRelativePositionsForSingleMatch )
{
    CigarMatch cigar( 1 );
    auto actualPositions = cigar.getRelativeRefPositions( 123 );
    BOOST_CHECK_EQUAL( actualPositions.size(), 1 );
    BOOST_CHECK_EQUAL( actualPositions.front(), 123 );
}

BOOST_AUTO_TEST_CASE( shouldGetTheCorrectRelativePositionsForThreeMatches )
{
    CigarMatch cigar( 3 );
    auto actualPositions = cigar.getRelativeRefPositions( 123 );
    std::vector< int64_t > expectedPositions{123, 124, 125};
    BOOST_CHECK_EQUAL_COLLECTIONS( actualPositions.begin(), actualPositions.end(), expectedPositions.begin(),
                                   expectedPositions.end() );
}

BOOST_AUTO_TEST_CASE( shouldGetTheCorrectRelativePositionsForThreeDeletions )
{
    CigarDeletion cigar( 3 );
    auto actualPositions = cigar.getRelativeRefPositions( 123 );
    std::vector< int64_t > expectedPositions{};
    BOOST_CHECK_EQUAL_COLLECTIONS( actualPositions.begin(), actualPositions.end(), expectedPositions.begin(),
                                   expectedPositions.end() );
}

BOOST_AUTO_TEST_CASE( shouldGetTheCorrectRelativePositionsForThreeInsertions )
{
    CigarInsertion cigar( 3 );
    auto actualPositions = cigar.getRelativeRefPositions( 123 );
    std::vector< int64_t > expectedPositions{emptyPos, emptyPos, emptyPos};
    BOOST_CHECK_EQUAL_COLLECTIONS( actualPositions.begin(), actualPositions.end(), expectedPositions.begin(),
                                   expectedPositions.end() );
}

BOOST_AUTO_TEST_CASE( shouldGetTheCorrectRelativePositionsForThreeSkips )
{
    CigarPadding cigar( 3 );
    auto actualPositions = cigar.getRelativeRefPositions( 123 );
    std::vector< int64_t > expectedPositions{};
    BOOST_CHECK_EQUAL_COLLECTIONS( actualPositions.begin(), actualPositions.end(), expectedPositions.begin(),
                                   expectedPositions.end() );
}

BOOST_AUTO_TEST_CASE( shouldGetTheCorrectRelativePositionsForThreePaddings )
{
    CigarPadding cigar( 3 );
    auto actualPositions = cigar.getRelativeRefPositions( 123 );
    std::vector< int64_t > expectedPositions{};
    BOOST_CHECK_EQUAL_COLLECTIONS( actualPositions.begin(), actualPositions.end(), expectedPositions.begin(),
                                   expectedPositions.end() );
}

BOOST_AUTO_TEST_CASE( shouldGetTheCorrectRelativePositionsForThreeSoftClips )
{
    CigarSoftClip cigar( 3 );
    auto actualPositions = cigar.getRelativeRefPositions( 123 );
    std::vector< int64_t > expectedPositions{emptyPos, emptyPos, emptyPos};
    BOOST_CHECK_EQUAL_COLLECTIONS( actualPositions.begin(), actualPositions.end(), expectedPositions.begin(),
                                   expectedPositions.end() );
}

BOOST_AUTO_TEST_CASE( shouldGetTheCorrectRelativePositionsForThreeHardClips )
{
    CigarHardClip cigar( 3 );
    auto actualPositions = cigar.getRelativeRefPositions( 123 );
    std::vector< int64_t > expectedPositions{};
    BOOST_CHECK_EQUAL_COLLECTIONS( actualPositions.begin(), actualPositions.end(), expectedPositions.begin(),
                                   expectedPositions.end() );
}

// Test should move offsets as follows.

BOOST_AUTO_TEST_CASE( shouldMoveOffsetsForMatch )
{
    auto offsets = std::make_shared< Offsets >( 0, 0 );
    CigarMatch cigar( 3 );
    cigar.moveOffsets( offsets );

    BOOST_CHECK_EQUAL( offsets->read, 3 );
    BOOST_CHECK_EQUAL( offsets->ref, 3 );
}

BOOST_AUTO_TEST_CASE( shouldMoveOffsetsForDeletion )
{
    auto offsets = std::make_shared< Offsets >( 0, 0 );
    CigarDeletion cigar( 3 );
    cigar.moveOffsets( offsets );

    BOOST_CHECK_EQUAL( offsets->read, 0 );
    BOOST_CHECK_EQUAL( offsets->ref, 3 );
}

BOOST_AUTO_TEST_CASE( shouldMoveOffsetsForInsertion )
{
    auto offsets = std::make_shared< Offsets >( 0, 0 );
    CigarInsertion cigar( 3 );
    cigar.moveOffsets( offsets );

    BOOST_CHECK_EQUAL( offsets->read, 3 );
    BOOST_CHECK_EQUAL( offsets->ref, 0 );
}

BOOST_AUTO_TEST_CASE( shouldMoveOffsetsForSkip )
{
    auto offsets = std::make_shared< Offsets >( 0, 0 );
    CigarSkip cigar( 3 );
    cigar.moveOffsets( offsets );

    BOOST_CHECK_EQUAL( offsets->read, 0 );
    BOOST_CHECK_EQUAL( offsets->ref, 3 );
}

BOOST_AUTO_TEST_CASE( shouldMoveOffsetsForPadding )
{
    auto offsets = std::make_shared< Offsets >( 0, 0 );
    CigarPadding cigar( 3 );
    cigar.moveOffsets( offsets );

    BOOST_CHECK_EQUAL( offsets->read, 0 );
    BOOST_CHECK_EQUAL( offsets->ref, 0 );
}

BOOST_AUTO_TEST_CASE( shouldMoveOffsetsForHardClip )
{
    auto offsets = std::make_shared< Offsets >( 0, 0 );
    CigarHardClip cigar( 3 );
    cigar.moveOffsets( offsets );

    BOOST_CHECK_EQUAL( offsets->read, 0 );
    BOOST_CHECK_EQUAL( offsets->ref, 0 );
}

BOOST_AUTO_TEST_CASE( shouldMoveOffsetsForSoftClip )
{
    auto offsets = std::make_shared< Offsets >( 0, 0 );
    CigarSoftClip cigar( 3 );
    cigar.moveOffsets( offsets );

    BOOST_CHECK_EQUAL( offsets->read, 3 );
    BOOST_CHECK_EQUAL( offsets->ref, 0 );
}

// Test get variants from cigar items.

BOOST_AUTO_TEST_CASE( shouldGetVariantsForMatch )
{
    const auto refSeq = std::make_shared< ReferenceSequence >( Region( "1", 0, 16 ), "ATATACTATATATATA" );
    const std::string readSeq = "AAATACTATATATATA";
    const std::string readQual = "QQQQQQQQQQQQQQQQ";

    auto variantGenerationData = std::make_shared< VariantGenerationData >( refSeq, 0, readSeq );
    auto offsets = std::make_shared< Offsets >( 0, 0 );
    CigarMatch cigar( 10 );

    auto variants = cigar.getVariants( variantGenerationData, offsets );
    BOOST_CHECK_EQUAL( variants.size(), 1 );
    BOOST_CHECK( checkVariantInVector( std::vector< varPtr_t >( variants.begin(), variants.end() ),
                                       std::make_shared< Variant >( refSeq, Region( "1", 1, 2 ), "A" ) ) );
}

BOOST_AUTO_TEST_CASE( shouldGetVariantsForDeletion )
{
    const Region region( "1", 0, 12 );
    const auto refSeq = std::make_shared< ReferenceSequence >(
        region, std::string( "G" ) + std::string( 10, 'A' ) + std::string( "G" ) );
    const std::string readSeq = "GG";
    const std::string readQual = "QQ";

    auto variantGenerationData = std::make_shared< VariantGenerationData >( refSeq, 0, readSeq );
    auto offsets = std::make_shared< Offsets >( 1, 1 );
    CigarDeletion cigar( 10 );

    auto variants = cigar.getVariants( variantGenerationData, offsets );
    BOOST_CHECK_EQUAL( variants.size(), 1 );
    BOOST_CHECK( checkVariantInVector( std::vector< varPtr_t >( variants.begin(), variants.end() ),
                                       std::make_shared< Variant >( refSeq, Region( "1", 1, 11 ), "" ) ) );
}

BOOST_AUTO_TEST_CASE( shouldNotGetVariantsForDeletionAtReadEnd )
{
    const std::string refSeqString = std::string( "G" ) + std::string( 1, 'A' );
    const Region region( "1", 0, refSeqString.size() );
    const auto refSeq = std::make_shared< ReferenceSequence >( region, refSeqString );
    const std::string readSeq = "G";
    const std::string readQual = "Q";

    auto variantGenerationData = std::make_shared< VariantGenerationData >( refSeq, 0, readSeq );
    auto offsets = std::make_shared< Offsets >( 1, 1 );
    CigarDeletion cigar( 1 );

    auto variants = cigar.getVariants( variantGenerationData, offsets );
    BOOST_CHECK_EQUAL( variants.size(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldNotGetVariantsForDeletionAtReadStart )
{
    const std::string refSeqString = std::string( 1, 'A' ) + std::string( "G" );
    const Region region( "1", 0, refSeqString.size() );
    const auto refSeq = std::make_shared< ReferenceSequence >( region, refSeqString );
    const std::string readSeq = "G";
    const std::string readQual = "Q";

    auto variantGenerationData = std::make_shared< VariantGenerationData >( refSeq, 0, readSeq );
    auto offsets = std::make_shared< Offsets >( 0, 0 );
    CigarDeletion cigar( 1 );

    auto variants = cigar.getVariants( variantGenerationData, offsets );
    BOOST_CHECK_EQUAL( variants.size(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldGetVariantsForInsertion )
{
    const std::string refSeqString = "GG";
    const Region region( "1", 0, refSeqString.size() );
    const auto refSeq = std::make_shared< ReferenceSequence >( region, refSeqString );
    const std::string readSeq = std::string( "G" ) + std::string( 10, 'A' ) + std::string( "G" );
    const std::string readQual = std::string( 12, 'A' );

    auto variantGenerationData = std::make_shared< VariantGenerationData >( refSeq, 0, readSeq );
    auto offsets = std::make_shared< Offsets >( 1, 1 );
    CigarInsertion cigar( 10 );

    auto variants = cigar.getVariants( variantGenerationData, offsets );
    BOOST_CHECK_EQUAL( variants.size(), 1 );
    BOOST_CHECK(
        checkVariantInVector( std::vector< varPtr_t >( variants.begin(), variants.end() ),
                              std::make_shared< Variant >( refSeq, Region( "1", 1, 1 ), std::string( 10, 'A' ) ) ) );
}

BOOST_AUTO_TEST_CASE( shouldGetVariantsForInsertionAtReadStart )
{
    const std::string refSeqString = "G";
    const Region region( "1", 0, refSeqString.size() );
    const auto refSeq = std::make_shared< ReferenceSequence >( region, refSeqString );
    const std::string readSeq = std::string( 1, 'A' ) + std::string( "G" );
    const std::string readQual = std::string( readSeq.size(), 'Q' );

    auto variantGenerationData = std::make_shared< VariantGenerationData >( refSeq, 0, readSeq );
    auto offsets = std::make_shared< Offsets >( 0, 0 );
    CigarInsertion cigar( 1 );

    auto variants = cigar.getVariants( variantGenerationData, offsets );
    BOOST_CHECK_EQUAL( variants.size(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldGetVariantsForInsertionAtReadEnd )
{
    const std::string refSeqString = "G";
    const Region region( "1", 0, refSeqString.size() );
    const auto refSeq = std::make_shared< ReferenceSequence >( region, refSeqString );
    const std::string readSeq = std::string( "G" ) + std::string( 1, 'A' );
    const std::string readQual = std::string( readSeq.size(), 'A' );

    auto variantGenerationData = std::make_shared< VariantGenerationData >( refSeq, 0, readSeq );
    auto offsets = std::make_shared< Offsets >( 1, 1 );
    CigarInsertion cigar( 1 );

    auto variants = cigar.getVariants( variantGenerationData, offsets );
    BOOST_CHECK_EQUAL( variants.size(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldGetEmptyVectorForHardClipCigar )
{
    CigarHardClip cigar( 10 );
    auto variants = cigar.getVariants( nullptr, nullptr );
    BOOST_CHECK_EQUAL( variants.size(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldGetEmptyVectorForSoftClipCigar )
{
    CigarSoftClip cigar( 10 );
    auto variants = cigar.getVariants( nullptr, nullptr );
    BOOST_CHECK_EQUAL( variants.size(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldGetEmptyVectorForPaddingCigar )
{
    CigarPadding cigar( 10 );
    auto variants = cigar.getVariants( nullptr, nullptr );
    BOOST_CHECK_EQUAL( variants.size(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldGetEmptyVectorForSkipCigar )
{
    CigarSkip cigar( 10 );
    auto variants = cigar.getVariants( nullptr, nullptr );
    BOOST_CHECK_EQUAL( variants.size(), 0 );
}
