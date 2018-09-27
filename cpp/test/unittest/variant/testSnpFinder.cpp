// All content Copyright (C) 2018 Genomics plc
#include "variant/variantGenerator.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "variant/type/variant.hpp"
#include "variant/snpFinder.hpp"
#include "unittest/vcf/VCFTestUtils.hpp"
#include "alignment/cigar.hpp"
#include "caller/region.hpp"

using echidna::variant::SNPFinder;
using echidna::variant::varPtr_t;
using echidna::variant::Variant;
using echidna::alignment::offsetsPtr_t;
using echidna::alignment::Offsets;
using echidna::variant::VariantGenerationData;
using echidna::alignment::Cigar;
using echidna::caller::Region;
using echidna::utils::ReferenceSequence;

//-------------------------------------------------------------------------------------------------

echidna::io::Read constructHighQualMatchReadWithSeqForSNPTesting( std::string seq, char baseQual = 60 )
{
    std::size_t seqLen = seq.length();
    std::string qual( seqLen, baseQual );
    size_t startPos = 0;
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( "1", startPos, startPos + seqLen + 20 ), std::string( seqLen + 20, 'A' ) );

    return echidna::io::Read( seq, qual, "testId", Cigar( std::to_string( seqLen ) + "M" ), 0, startPos, 0, 0, 0, 0, 0,
                              referenceSequence );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( shouldFindSNPsAtEdges )
{
    std::string refSeqString( "AAAAA" );
    auto read = constructHighQualMatchReadWithSeqForSNPTesting( "TAAAT" );

    Region region( "chr1", read.getStartPos(), read.getAlignedEndPos() );
    const auto refSeq = std::make_shared< ReferenceSequence >( region, refSeqString );

    auto variantGenerationData = std::make_shared< VariantGenerationData >( refSeq, read );
    SNPFinder finder( variantGenerationData );

    auto variants =
        finder.findSNPsInReadSegment( std::make_shared< Offsets >( 0, 0 ), read.getEndPos() - read.getStartPos() );

    std::vector< varPtr_t > vecVariants( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( variants.size(), 2 );
    BOOST_CHECK(
        checkVariantInVector( vecVariants, std::make_shared< Variant >( refSeq, Region( "chr1", 0, 1 ), "T" ) ) );
    BOOST_CHECK(
        checkVariantInVector( vecVariants, std::make_shared< Variant >( refSeq, Region( "chr1", 4, 5 ), "T" ) ) );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( shouldCallSNPsWithReadOffsetAndZeroMargin )
{
    std::string refSeqString( "ABCDE" );
    auto read = constructHighQualMatchReadWithSeqForSNPTesting( "GGCTGCDE" );

    auto offsets = std::make_shared< Offsets >( 3, 0 );  // read, ref
    Region region( "chr1", read.getStartPos(), read.getStartPos() + refSeqString.size() );
    const auto refSeq = std::make_shared< ReferenceSequence >( region, refSeqString );

    auto variantGenerationData = std::make_shared< VariantGenerationData >( refSeq, read );
    SNPFinder finder( variantGenerationData );

    auto variants = finder.findSNPsInReadSegment( offsets, read.getEndPos() - read.getStartPos() );

    std::vector< varPtr_t > vecVariants( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( variants.size(), 2 );
    BOOST_CHECK(
        checkVariantInVector( vecVariants, std::make_shared< Variant >( refSeq, Region( "chr1", 0, 1 ), "T" ) ) );
    BOOST_CHECK(
        checkVariantInVector( vecVariants, std::make_shared< Variant >( refSeq, Region( "chr1", 1, 2 ), "G" ) ) );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( shouldCallSNPsWithRefOffsetAndZeroMargin )
{
    std::string refSeqString( "ABCDE" );
    auto read = constructHighQualMatchReadWithSeqForSNPTesting( "GCGT" );

    auto offsets = std::make_shared< Offsets >( 0, 3 );  // read, ref
    Region region( "chr1", read.getStartPos(), read.getStartPos() + refSeqString.size() );
    const auto refSeq = std::make_shared< ReferenceSequence >( region, refSeqString );

    auto variantGenerationData = std::make_shared< VariantGenerationData >( refSeq, read );
    SNPFinder finder( variantGenerationData );

    auto variants = finder.findSNPsInReadSegment( offsets, read.getEndPos() - read.getStartPos() );

    std::vector< varPtr_t > vecVariants( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( variants.size(), 2 );
    BOOST_CHECK(
        checkVariantInVector( vecVariants, std::make_shared< Variant >( refSeq, Region( "chr1", 3, 4 ), "G" ) ) );
    BOOST_CHECK(
        checkVariantInVector( vecVariants, std::make_shared< Variant >( refSeq, Region( "chr1", 4, 5 ), "C" ) ) );
}
//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( shouldCallSNPsWithLength )
{
    std::string refSeqString( "AAAAA" );
    auto read = constructHighQualMatchReadWithSeqForSNPTesting( "GCGTG" );

    auto offsets = std::make_shared< Offsets >( 0, 0 );  // read, ref
    size_t length = 2;
    Region region( "chr1", read.getStartPos(), read.getStartPos() + refSeqString.size() );
    const auto refSeq = std::make_shared< ReferenceSequence >( region, refSeqString );

    auto variantGenerationData = std::make_shared< VariantGenerationData >( refSeq, read );
    SNPFinder finder( variantGenerationData );

    auto variants = finder.findSNPsInReadSegment( offsets, length );

    std::vector< varPtr_t > vecVariants( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( variants.size(), 2 );
    BOOST_CHECK(
        checkVariantInVector( vecVariants, std::make_shared< Variant >( refSeq, Region( "chr1", 0, 1 ), "G" ) ) );
    BOOST_CHECK(
        checkVariantInVector( vecVariants, std::make_shared< Variant >( refSeq, Region( "chr1", 1, 2 ), "C" ) ) );
}
