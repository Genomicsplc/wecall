#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "io/read.hpp"
#include "io/readRange.hpp"
#include "io/readDataSet.hpp"
#include "alignment/cigar.hpp"
#include "alignment/cigarItems.hpp"
#include "common.hpp"

using echidna::io::Read;
using echidna::alignment::Cigar;
using echidna::utils::ReferenceSequence;
using echidna::utils::BasePairSequence;
using echidna::utils::QualitySequence;
using echidna::caller::Region;
using echidna::variant::Variant;
using echidna::variant::Breakpoint;
using echidna::io::ReadDataset;

std::shared_ptr< Read > constructHighQualMatchReadWithLength( size_t seqLength, int64_t startPos )
{
    std::string seq( seqLength, 'T' );
    std::string qual( seqLength, 60 );

    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >(
        Region( "1", startPos - 1, startPos + seqLength + 1 ), std::string( seqLength + 2, 'A' ) );
    return std::make_shared< Read >( seq, qual, "testId", Cigar( std::to_string( seqLength ) + "M" ), 0, startPos, 0, 0,
                                     0, 0, 0, refSequence );
}

BOOST_AUTO_TEST_CASE( testReadsStartEndComputationForEmptyList )
{
    ReadDataset readDataset( {""}, Region( "1", 0, 10 ) );

    const auto startEnd = echidna::io::readsAlignedStartEnd( readDataset.getAllReads( 0 ).at( "" ) );
    BOOST_CHECK_EQUAL( startEnd.first, echidna::alignment::noPos );
    BOOST_CHECK_EQUAL( startEnd.second, echidna::alignment::noPos );
}

BOOST_AUTO_TEST_CASE( testReadsStartEndComputation )
{
    ReadDataset readDataset( {""}, Region( "1", 0, 10 ) );
    auto refSequenceLen = 100;
    auto refSequence = std::make_shared< ReferenceSequence >( Region( "1", -10, refSequenceLen - 10 ),
                                                              std::string( refSequenceLen, 'A' ) );
    readDataset.insertRead( "", std::make_shared< Read >( std::string( 8, 'A' ), std::string( 8, 'Q' ), "1",
                                                          Cigar( "8M" ), 0, 3, 0, 0, 0, 0, 0, refSequence ) );
    readDataset.insertRead( "", std::make_shared< Read >( std::string( 11, 'A' ), std::string( 11, 'Q' ), "2",
                                                          Cigar( "11M" ), 0, 0, 0, 0, 0, 0, 0, refSequence ) );
    readDataset.insertRead( "", std::make_shared< Read >( std::string( 6, 'A' ), std::string( 6, 'Q' ), "3",
                                                          Cigar( "6M" ), 0, 4, 0, 0, 0, 0, 0, refSequence ) );

    const auto startEnd = echidna::io::readsAlignedStartEnd( readDataset.getAllReads( 0 ).at( "" ) );
    BOOST_CHECK_EQUAL( startEnd.first, 0 );
    BOOST_CHECK_EQUAL( startEnd.second, 11 );
}

BOOST_AUTO_TEST_CASE( testMaxAlignedReadsLengthComputationEmptyTrees )
{
    ReadDataset readDataset( {"A", "B", "C", "D"}, Region( "1", 0, 10 ) );
    const auto reads = readDataset.getAllReads( 0 );

    BOOST_CHECK_EQUAL( echidna::io::perSampleMaxAlignedReadLength( reads ), 0 );
    BOOST_CHECK_EQUAL( echidna::io::perSampleMaxReadCigarLength( reads ), 0 );
}

BOOST_AUTO_TEST_CASE( testMaxAlignedReadsLengthComputation )
{
    ReadDataset readDataset( {"A", "B", "C", "D"}, Region( "1", 0, 100 ) );
    readDataset.insertRead( "A", constructHighQualMatchReadWithLength( 8, 0 ) );
    readDataset.insertRead( "B", constructHighQualMatchReadWithLength( 11, 0 ) );
    readDataset.insertRead( "B", constructHighQualMatchReadWithLength( 6, 4 ) );
    readDataset.insertRead( "C", constructHighQualMatchReadWithLength( 80, 20 ) );
    readDataset.insertRead( "C", constructHighQualMatchReadWithLength( 82, 17 ) );

    const auto reads = readDataset.getAllReads( 0 );

    BOOST_CHECK_EQUAL( echidna::io::perSampleMaxAlignedReadLength( reads ), 99 - 17 );
    BOOST_CHECK_EQUAL( echidna::io::perSampleMaxReadCigarLength( reads ), 99 - 17 );
}

BOOST_AUTO_TEST_CASE( testMaxReadsLengthComputationEmptyTrees )
{
    ReadDataset readDataset( {"A", "B", "C", "D"}, Region( "1", 0, 10 ) );
    const auto reads = readDataset.getAllReads( 0 );

    BOOST_CHECK_EQUAL( echidna::io::perSampleMaxReadLength( reads ), 0 );
    BOOST_CHECK_EQUAL( echidna::io::perSampleMaxReadCigarLength( reads ), 0 );
}

BOOST_AUTO_TEST_CASE( testMaxReadsLengthComputation )
{
    ReadDataset readDataset( {"A", "B", "C", "D"}, Region( "1", 0, 100 ) );
    readDataset.insertRead( "A", constructHighQualMatchReadWithLength( 8, 0 ) );
    readDataset.insertRead( "B", constructHighQualMatchReadWithLength( 11, 0 ) );
    readDataset.insertRead( "B", constructHighQualMatchReadWithLength( 6, 4 ) );
    readDataset.insertRead( "C", constructHighQualMatchReadWithLength( 80, 20 ) );
    readDataset.insertRead( "C", constructHighQualMatchReadWithLength( 82, 17 ) );

    const auto reads = readDataset.getAllReads( 0 );

    BOOST_CHECK_EQUAL( echidna::io::perSampleMaxReadLength( reads ), 99 - 17 );
    BOOST_CHECK_EQUAL( echidna::io::perSampleMaxReadCigarLength( reads ), 99 - 17 );
}

BOOST_AUTO_TEST_CASE( testGetReferencePositionsSimpleCase )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );
    Read read( "WHY", "NOT", "test", Cigar( "3M" ), 0, 0, BAM_FPROPER_PAIR, 100, 200, 0, 200, refSequence );

    auto positions = read.getReferencePositions();
    auto expected_positions = {0, 1, 2};
    BOOST_CHECK_EQUAL_COLLECTIONS( positions.begin(), positions.end(), expected_positions.begin(),
                                   expected_positions.end() );
}

BOOST_AUTO_TEST_CASE( testGetReferencePositionsFromRead )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );
    Read read( "WHY", "NOT", "test", Cigar( "1M10D1I1M" ), 0, 0, BAM_FPROPER_PAIR, 100, 200, 0, 200, refSequence );

    auto positions = read.getReferencePositions();
    auto expected_positions = {0, -1, 11};
    BOOST_CHECK_EQUAL_COLLECTIONS( positions.begin(), positions.end(), expected_positions.begin(),
                                   expected_positions.end() );
}

BOOST_AUTO_TEST_CASE( shouldtrimReadOfShortFragmentAtEnd )
{
    std::size_t length = 10;
    std::string seq( length, 'A' );
    std::string qual( length, 'Q' );
    std::string strCig = "10M";
    int64_t startPos = 0;
    int64_t flag = BAM_FPROPER_PAIR;

    int64_t insertSize = 9;

    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", -10, 10 ), std::string( 20, 'A' ) );

    Read read( seq, qual, "Fwd", Cigar( strCig ), 0, startPos, flag, 100, insertSize, 0, 200, refSequence );
    read.trimReadOfShortFragment();

    auto lengthToTrim = length - int64_to_sizet( insertSize );

    std::string expectedQualities =
        qual.substr( 0, qual.size() - lengthToTrim ) + std::string( lengthToTrim, constants::minAllowedQualityScore );
    BOOST_CHECK_EQUAL( read.getQualities(), expectedQualities );
}

BOOST_AUTO_TEST_CASE( shouldNotTrimReadIfInsertSizeIsGreaterThanReadLength )
{
    std::size_t length = 10;
    std::string seq( length, 'A' );
    std::string qual( length, 'Q' );
    std::string strCig = "10M";
    int64_t startPos = 0;
    int64_t flag = BAM_FPROPER_PAIR + BAM_FREVERSE;
    int64_t insertSize = length + 100000000;
    int64_t mappingQuality = 100;
    int64_t mateStartPos = startPos + insertSize - length;

    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );
    Read read( seq, qual, "Rev", Cigar( strCig ), 0, startPos, flag, mappingQuality, insertSize, 0, mateStartPos,
               refSequence );
    read.trimReadOfShortFragment();

    std::string expectedQualities = qual;
    BOOST_CHECK_EQUAL( read.getQualities(), expectedQualities );
}

BOOST_AUTO_TEST_CASE( shouldtrimReadOfShortFragmentReverseReadAtStart )
{
    std::size_t length = 10;
    std::string seq( length, 'A' );
    std::string qual( length, 'Q' );
    std::string strCig = "10M";
    int64_t mappingQuality = 100;
    int64_t startPos = 0;
    int64_t insertSize = 9;
    int64_t mateStartPos = startPos + insertSize - length;
    int64_t flag = BAM_FPROPER_PAIR + BAM_FREVERSE;

    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );
    Read read( seq, qual, "Rev", Cigar( strCig ), 0, startPos, flag, mappingQuality, insertSize, 0, mateStartPos,
               refSequence );
    read.trimReadOfShortFragment();

    auto lengthToTrim = length - int64_to_sizet( insertSize );
    std::string expectedQualities = std::string( lengthToTrim, constants::minAllowedQualityScore ) +
                                    qual.substr( lengthToTrim, qual.size() - lengthToTrim );
    BOOST_CHECK_EQUAL( read.getQualities(), expectedQualities );
}

BOOST_AUTO_TEST_CASE( shouldTrimOverlapOfForwardRead2 )
{
    std::size_t length = 10;
    std::string seq( length, 'A' );
    std::string qual( length, 'Q' );
    std::string strCig = "10M";
    int64_t startPos = 0;
    int64_t flag = BAM_FPROPER_PAIR;

    std::size_t overlapLength = 1;
    int64_t insertSize = 2 * static_cast< int64_t >( length ) - static_cast< int64_t >( overlapLength );

    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );
    Read read( seq, qual, "Fwd", Cigar( strCig ), 0, startPos, flag, 100, insertSize, 0, 200, refSequence );
    read.trimOverlap();

    std::string expectedQualities =
        qual.substr( 0, qual.size() - overlapLength ) + std::string( overlapLength, constants::minAllowedQualityScore );

    BOOST_CHECK_EQUAL( read.getQualities(), expectedQualities );
}

BOOST_AUTO_TEST_CASE( shouldNotTrimOverlapOfForwardRead1 )
{
    std::size_t length = 10;
    std::string seq( length, 'A' );
    std::string qual( length, 'Q' );
    std::string strCig = "10M";
    int64_t startPos = 0;
    int64_t flag = BAM_FPROPER_PAIR + BAM_FREAD1;

    std::size_t overlapLength = 1;
    int64_t insertSize = 2 * static_cast< int64_t >( length ) - static_cast< int64_t >( overlapLength );

    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );
    Read read( seq, qual, "Fwd", Cigar( strCig ), 0, startPos, flag, 100, insertSize, 0, 200, refSequence );
    read.trimOverlap();

    std::string expectedQualities = qual;

    BOOST_CHECK_EQUAL( read.getQualities(), expectedQualities );
}

BOOST_AUTO_TEST_CASE( shouldTrimOverlapOfReverseRead2 )
{
    std::size_t length = 10;
    std::string seq( length, 'A' );
    std::string qual( length, 'Q' );
    std::string strCig = "10M";
    int64_t startPos = 0;
    int64_t flag = BAM_FPROPER_PAIR + BAM_FREVERSE;

    std::size_t overlapLength = 1;
    int64_t insertSize = 2 * static_cast< int64_t >( length ) - static_cast< int64_t >( overlapLength );

    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );
    Read read( seq, qual, "Rev", Cigar( strCig ), 0, startPos, flag, 100, insertSize, 0, 200, refSequence );
    read.trimOverlap();

    std::string expectedQualities = std::string( overlapLength, constants::minAllowedQualityScore ) +
                                    qual.substr( overlapLength, qual.size() - overlapLength );

    BOOST_CHECK_EQUAL( read.getQualities(), expectedQualities );
}

BOOST_AUTO_TEST_CASE( shouldGetWholeReadSpanIfIntervalMatchesReadSpanInRef )
{
    Cigar cigar( "2M2I2M2D2M" );
    const auto startPos = 100L;
    const std::string seq( 8, 'A' );
    const std::string qual( 8, 'Q' );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    BOOST_CHECK_EQUAL( testRead.getAlignedEndPos(), 108L );

    echidna::utils::Interval inputRefInterval( startPos, 108L );
    echidna::utils::Interval expectedResult( 0L, testRead.getLength() );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldGetWholeReadSpanIfIntervalContainsReadSpanInRef )
{
    Cigar cigar( "2M2I2M2D2M" );
    const auto startPos = 100L;
    const std::string seq( 8, 'A' );
    const std::string qual( 8, 'Q' );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::utils::Interval inputRefInterval( startPos - 10, testRead.getAlignedEndPos() + 10 );
    echidna::utils::Interval expectedResult( 0L, testRead.getLength() );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldGetStartOfReadSpanIfIntervalPreceedsReadSpanInRef )
{
    Cigar cigar( "2M2I2M2D2M" );
    const auto startPos = 100L;
    const std::string seq( 8, 'A' );
    const std::string qual( 8, 'Q' );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::utils::Interval inputRefInterval( startPos - 10, startPos );
    echidna::utils::Interval expectedResult( 0L, 0L );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldGetEndOfReadSpanIfIntervalFollowsReadSpanInRef )
{
    Cigar cigar( "2M2I2M2D2M" );
    const auto startPos = 100L;
    const std::string seq( 8, 'A' );
    const std::string qual( 8, 'Q' );
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 100, 110 ), std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::utils::Interval inputRefInterval( testRead.getAlignedEndPos(), testRead.getAlignedEndPos() + 10 );
    echidna::utils::Interval expectedResult( testRead.getLength(), testRead.getLength() );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldGetEmptyIntervalIfInputIntervalIsEmptyAndReadFlatAligned )
{
    Cigar cigar( "2M2I2M2D2M" );
    const auto startPos = 100L;
    const std::string seq( 8, 'A' );
    const std::string qual( 8, 'Q' );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::utils::Interval inputRefInterval( startPos + 1L, startPos + 1L );
    echidna::utils::Interval expectedResult( 1L, 1L );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldGetIntervalCorrespondingToInsertionInReadIfInputIntervalIsEmptyAtStartPos )
{
    Cigar cigar( "2M2I2M2D2M" );
    const auto startPos = 100L;
    const std::string seq( 8, 'A' );
    const std::string qual( 8, 'Q' );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::utils::Interval inputRefInterval( startPos + 2L, startPos + 2L );
    echidna::utils::Interval expectedResult( 2L, 4L );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldGetMatchingIntervalIfInputIntervalIsFlatAligned )
{
    Cigar cigar( "2M2I2M2D2M" );
    const auto startPos = 100L;
    const std::string seq( 8, 'A' );
    const std::string qual( 8, 'Q' );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::utils::Interval inputRefInterval( startPos + 3L, startPos + 4L );
    echidna::utils::Interval expectedResult( 5L, 6L );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldGetEmptyIntervalForEmptyRefIntervalCorrespondingToDeletion )
{
    Cigar cigar( "2M2I2M2D2M" );
    const auto startPos = 100L;
    const std::string seq( 8, 'A' );
    const std::string qual( 8, 'Q' );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::utils::Interval inputRefInterval( startPos + 5L, startPos + 5L );
    echidna::utils::Interval expectedResult( 6L, 6L );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldGetEmptyIntervalForNonEmptyRefIntervalCorrespondingToDeletion )
{
    Cigar cigar( "2M2I2M2D2M" );
    const auto startPos = 100L;
    const std::string seq( 8, 'A' );
    const std::string qual( 8, 'Q' );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::utils::Interval inputRefInterval( startPos + 4L, startPos + 6L );
    echidna::utils::Interval expectedResult( 6L, 6L );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

// indels at edges of reads

BOOST_AUTO_TEST_CASE( shouldGetInsertionAtStartOfReadForNonOverlappingInterval )
{
    Cigar cigar( "4I4M" );
    const auto startPos = 100L;
    const std::string seq( 8, 'A' );
    const std::string qual( 8, 'Q' );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::utils::Interval inputRefInterval( startPos + 0L, startPos + 0L );
    echidna::utils::Interval expectedResult( 0L, 4L );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldIgnoreDeletionAtStartOfReadForNonOverlappingInterval )
{
    Cigar cigar( "4D4M" );
    const auto startPos = 100L;
    const std::string seq( 4, 'A' );
    const std::string qual( 4, 'Q' );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::utils::Interval inputRefInterval( startPos + 0L, startPos + 0L );
    echidna::utils::Interval expectedResult( 0L, 0L );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldGetInsertionAtEndOfReadForNonOverlappingInterval )
{
    Cigar cigar( "4M4I" );
    const auto startPos = 100L;
    const std::string seq( 8, 'A' );
    const std::string qual( 8, 'Q' );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::utils::Interval inputRefInterval( testRead.getAlignedEndPos(), testRead.getAlignedEndPos() );
    echidna::utils::Interval expectedResult( 4L, 8L );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldIgnoreDeletionAtEndOfReadForNonOverlappingInterval )
{
    Cigar cigar( "4M4D" );
    const auto startPos = 100L;
    const std::string seq( 4, 'A' );
    const std::string qual( 4, 'Q' );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::utils::Interval inputRefInterval( testRead.getAlignedEndPos(), testRead.getAlignedEndPos() );
    echidna::utils::Interval expectedResult( 4L, 4L );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

// intervals adjacent to indels

BOOST_AUTO_TEST_CASE( shouldGetInsertionOnLeftOfQueryInterval )
{
    Cigar cigar( "1M4I3M" );
    const auto startPos = 100L;
    const std::string seq( 8, 'A' );
    const std::string qual( 8, 'Q' );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::utils::Interval inputRefInterval( startPos + 1L, startPos + 2L );
    echidna::utils::Interval expectedResult( 1L, 6L );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldSkipDeletionOnLeftOfQueryInterval )
{
    Cigar cigar( "1M4D3M" );
    const auto startPos = 100L;
    const std::string seq( 4, 'A' );
    const std::string qual( 4, 'Q' );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::utils::Interval inputRefInterval( startPos + 5L, startPos + 6L );
    echidna::utils::Interval expectedResult( 1L, 2L );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldGetInsertionOnRightOfQueryInterval )
{
    Cigar cigar( "3M4I1M" );
    const auto startPos = 100L;
    const std::string seq( 8, 'A' );
    const std::string qual( 8, 'Q' );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::utils::Interval inputRefInterval( startPos + 2L, startPos + 3L );
    echidna::utils::Interval expectedResult( 2L, 7L );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldSkipDeletionOnRightOfQueryInterval )
{
    Cigar cigar( "3M4D1M" );
    const auto startPos = 100L;
    const std::string seq( 4, 'A' );
    const std::string qual( 4, 'Q' );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read testRead( seq, qual, "", cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::utils::Interval inputRefInterval( startPos + 2L, startPos + 3L );
    echidna::utils::Interval expectedResult( 2L, 3L );

    BOOST_CHECK_EQUAL( expectedResult, testRead.getIntervalInRead( inputRefInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldRetreiveSNPsFromRead )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "AAAAA" );
    const Read read( BasePairSequence( "TACG" ), QualitySequence( "QQQQ" ), "", Cigar( "4M" ), 0, 1, 0, 0, 0, 0, 0,
                     referenceSequence );
    const auto variants = read.getVariants();
    BOOST_REQUIRE_EQUAL( variants.size(), 3 );
    BOOST_CHECK_EQUAL( *variants[0], Variant( referenceSequence, Region( "1", 1, 2 ), "T" ) );
    BOOST_CHECK_EQUAL( *variants[1], Variant( referenceSequence, Region( "1", 3, 4 ), "C" ) );
    BOOST_CHECK_EQUAL( *variants[2], Variant( referenceSequence, Region( "1", 4, 5 ), "G" ) );
}

BOOST_AUTO_TEST_CASE( shouldRetrieveDeletionFromReadStart )
{
    bool skip = true;
    if ( skip )
    {
        return;
    }
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    const auto startPos = 1;
    const Read read( BasePairSequence( "A" ), QualitySequence( {10} ), "", Cigar( "4D1M" ), 0, startPos, 0, 0, 0, 0, 0,
                     referenceSequence );
    const auto variants = read.getVariants();
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    BOOST_CHECK_EQUAL( *variants[0], Variant( referenceSequence, Region( "1", 1, 5 ), "" ) );
}

BOOST_AUTO_TEST_CASE( shouldRetrieveDeletionFromReadEnd )
{
    bool skip = true;
    if ( skip )
    {
        return;
    }
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    const auto startPos = 1;
    const Read read( BasePairSequence( "A" ), QualitySequence( {10} ), "", Cigar( "1M4D" ), 0, startPos, 0, 0, 0, 0, 0,
                     referenceSequence );
    const auto variants = read.getVariants();
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    BOOST_CHECK_EQUAL( *variants[0], Variant( referenceSequence, Region( "1", 2, 6 ), "" ) );
}

BOOST_AUTO_TEST_CASE( shouldRetrieveDeletionFromReadMiddle )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    const auto start = 1;
    const Read read( BasePairSequence( "AA" ), QualitySequence( {10, 10} ), "", Cigar( "1M4D1M" ), 0, start, 0, 0, 0, 0,
                     0, referenceSequence );
    const auto variants = read.getVariants();
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    BOOST_CHECK_EQUAL( *variants[0], Variant( referenceSequence, Region( "1", 2, 6 ), "" ) );
}

BOOST_AUTO_TEST_CASE( shouldRetrieveInsertionFromReadStart )
{
    bool skip = true;
    if ( skip )
    {
        return;
    }
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    const auto startPos = 1;
    const Read read( BasePairSequence( "AAAAA" ), QualitySequence( 5, 10 ), "", Cigar( "4I1M" ), 0, startPos, 0, 0, 0,
                     0, 0, referenceSequence );
    const auto variants = read.getVariants();
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    BOOST_CHECK_EQUAL( *variants[0], Variant( referenceSequence, Region( "1", 1, 1 ), "AAAA" ) );
}

BOOST_AUTO_TEST_CASE( shouldRetrieveInsertionFromReadEnd )
{
    bool skip = true;
    if ( skip )
    {
        return;
    }
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    const auto startPos = 1;
    const Read read( BasePairSequence( "AAAAA" ), QualitySequence( 5, 10 ), "", Cigar( "1M4I" ), 0, startPos, 0, 0, 0,
                     0, 0, referenceSequence );
    const auto variants = read.getVariants();
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    BOOST_CHECK_EQUAL( *variants[0], Variant( referenceSequence, Region( "1", 2, 2 ), "AAAA" ) );
}

BOOST_AUTO_TEST_CASE( needsTwoCigarItemsToGetBreakpoints )
{
    const auto startPos = 1;
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read read( BasePairSequence( "AAAAA" ), QualitySequence( 5, 10 ), "", Cigar( "5S" ), 0, startPos, 0, 0, 0, 0,
                     0, refSequence );
    const auto breakpoints = read.getBreakpoints();
    BOOST_REQUIRE_EQUAL( breakpoints.size(), 0 );
}

BOOST_AUTO_TEST_CASE( shouldRetrieveBreakpointFromReadStart )
{
    const auto startPos = 1;
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read read( BasePairSequence( "AAAAA" ), QualitySequence( 5, 10 ), "", Cigar( "4S1M" ), 0, startPos, 0, 0, 0,
                     0, 0, refSequence );
    const auto breakpoints = read.getBreakpoints();
    BOOST_REQUIRE_EQUAL( breakpoints.size(), 1 );
    BOOST_CHECK_EQUAL( *breakpoints[0], Breakpoint( "1", startPos, false, "AAAA" ) );
}

BOOST_AUTO_TEST_CASE( shouldRetrieveBreakpointFromReadEnd )
{
    const auto startPos = 1;
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read read( BasePairSequence( "AAAAA" ), QualitySequence( 5, 10 ), "", Cigar( "1M4S" ), 0, startPos, 0, 0, 0,
                     0, 0, refSequence );
    const auto breakpoints = read.getBreakpoints();
    BOOST_REQUIRE_EQUAL( breakpoints.size(), 1 );
    BOOST_CHECK_EQUAL( *breakpoints[0], Breakpoint( "1", startPos + 1, true, "AAAA" ) );
}

BOOST_AUTO_TEST_CASE( testStoringOfMateRegion )
{
    const int32_t tid = 0;
    const int64_t startPos = 0;
    const int64_t flag = BAM_FPROPER_PAIR;
    int64_t mateStartPos = 200;

    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const Read read( BasePairSequence( "AAAAA" ), QualitySequence( 5, 10 ), "", Cigar( "1M4S" ), tid, startPos, flag, 0,
                     0, tid, mateStartPos, refSequence );
    const auto mateRegions = read.getBreakpoints()[0]->mateRegions();

    // Currently padding by 150 to left + using length of read as proxy for length of mate.
    std::set< Region > expected = {Region( "1", read.getMateIntervalInRef() )};

    BOOST_CHECK_EQUAL_COLLECTIONS( mateRegions.cbegin(), mateRegions.cend(), expected.cbegin(), expected.cend() );
}

BOOST_AUTO_TEST_CASE( testStoringOfMateRegionGetsNoubtIfMateUnmapped )
{
    const int32_t tid = 0;
    const int64_t startPos = 0;
    const int64_t flag = BAM_FPROPER_PAIR + BAM_FMUNMAP;
    int64_t mateStartPos = 200;

    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", startPos - 1, startPos + 9 ),
                                                                              std::string( 10, 'A' ) );
    const auto mateRegions = Read( BasePairSequence( "AAAAA" ), QualitySequence( 5, 10 ), "", Cigar( "1M4S" ), tid,
                                   startPos, flag, 0, 0, tid, mateStartPos, refSequence )
                                 .getBreakpoints()[0]
                                 ->mateRegions();

    std::set< Region > expected = {};

    BOOST_CHECK_EQUAL_COLLECTIONS( mateRegions.cbegin(), mateRegions.cend(), expected.cbegin(), expected.cend() );
}

BOOST_AUTO_TEST_CASE( shouldRetrieveInsertionFromReadMiddle )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    const auto start = 1;
    const Read read( BasePairSequence( "AAAAAA" ), QualitySequence( 6, 10 ), "", Cigar( "1M4I1M" ), 0, start, 0, 0, 0,
                     0, 0, referenceSequence );
    const auto variants = read.getVariants();
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    BOOST_CHECK_EQUAL( *variants[0], Variant( referenceSequence, Region( "1", 2, 2 ), "AAAA" ) );
}

BOOST_AUTO_TEST_CASE( shouldRetrieveSNPsEitherSideOfInsertion )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    const auto start = 1;
    const Read read( BasePairSequence( "TAAAAC" ), QualitySequence( 6, 10 ), "", Cigar( "1M4I1M" ), 0, start, 0, 0, 0,
                     0, 0, referenceSequence );
    const auto variants = read.getVariants();
    BOOST_REQUIRE_EQUAL( variants.size(), 3 );
    BOOST_CHECK_EQUAL( *variants[0], Variant( referenceSequence, Region( "1", 1, 2 ), "T" ) );
    BOOST_CHECK_EQUAL( *variants[1], Variant( referenceSequence, Region( "1", 2, 2 ), "AAAA" ) );
    BOOST_CHECK_EQUAL( *variants[2], Variant( referenceSequence, Region( "1", 2, 3 ), "C" ) );
}
