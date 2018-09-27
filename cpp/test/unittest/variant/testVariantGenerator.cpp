// All content Copyright (C) 2018 Genomics plc
#include "variant/variantGenerator.hpp"
#include "variant/type/variant.hpp"
#include "io/read.hpp"
#include "alignment/cigar.hpp"
#include "io/readDataSet.hpp"

#include "unittest/vcf/VCFTestUtils.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <string>

using namespace echidna::variant;
using namespace echidna::alignment;
using namespace echidna::io;
using ReadContainer = echidna::io::readIntervalTree_t;
using ReadIterator = ReadContainer::iterator;
using echidna::caller::Region;
using echidna::utils::ReferenceSequence;
using echidna::io::ReadDataset;

// TODO(ES): Breakpoints

BOOST_AUTO_TEST_CASE( testShouldRetrieveCorrectQualitiesForSNPs )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "CATGA" );

    auto var1 = std::make_shared< Variant >( alignedSequence, Region( "1", 0, 1 ), "G" );
    auto var2 = std::make_shared< Variant >( alignedSequence, Region( "1", 2, 3 ), "C" );
    auto var3 = std::make_shared< Variant >( alignedSequence, Region( "1", 4, 5 ), "T" );

    const auto qualities = echidna::variant::getVariantsReadBaseQualities( 0, "STEFT", {var1, var2, var3}, {} );

    const auto expectedResults = {
        static_cast< phred_t >( 'S' ), static_cast< phred_t >( 'E' ), static_cast< phred_t >( 'T' )};
    BOOST_CHECK_EQUAL_COLLECTIONS( expectedResults.begin(), expectedResults.end(), qualities.begin(), qualities.end() );
}

BOOST_AUTO_TEST_CASE( testShouldRetrieveCorrectQualitiesForSNPsAndDeletion )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "CATGA" );

    auto var1 = std::make_shared< Variant >( alignedSequence, Region( "1", 0, 1 ), "G" );
    auto var2 = std::make_shared< Variant >( alignedSequence, Region( "1", 2, 3 ), "" );
    auto var3 = std::make_shared< Variant >( alignedSequence, Region( "1", 4, 5 ), "T" );

    const auto qualities = echidna::variant::getVariantsReadBaseQualities( 0, "STEF", {var1, var2, var3}, {} );

    const auto expectedResults = {static_cast< phred_t >( 'S' ), 1000.0, static_cast< phred_t >( 'F' )};
    BOOST_CHECK_EQUAL_COLLECTIONS( expectedResults.begin(), expectedResults.end(), qualities.begin(), qualities.end() );
}

BOOST_AUTO_TEST_CASE( testShouldRetrieveCorrectQualitiesForSNPsAndInsertion )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "CATGA" );

    auto var1 = std::make_shared< Variant >( alignedSequence, Region( "1", 0, 1 ), "G" );
    auto var2 = std::make_shared< Variant >( alignedSequence, Region( "1", 2, 2 ), "C" );
    auto var3 = std::make_shared< Variant >( alignedSequence, Region( "1", 4, 5 ), "T" );

    const auto qualities = echidna::variant::getVariantsReadBaseQualities( 0, "STEFTH", {var1, var2, var3}, {} );

    const auto expectedResults = {
        static_cast< phred_t >( 'S' ), static_cast< phred_t >( 'E' ), static_cast< phred_t >( 'H' )};
    BOOST_CHECK_EQUAL_COLLECTIONS( expectedResults.begin(), expectedResults.end(), qualities.begin(), qualities.end() );
}

BOOST_AUTO_TEST_CASE( testShouldRetrieveCorrectQualitiesForDeletion )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "CATGA" );

    auto var1 = std::make_shared< Variant >( alignedSequence, Region( "1", 1, 4 ), "" );
    // Check alignment correctly altered
    auto var2 = std::make_shared< Variant >( alignedSequence, Region( "1", 4, 5 ), "T" );

    const auto qualities = echidna::variant::getVariantsReadBaseQualities( 0, "ST", {var1, var2}, {} );

    const auto expectedResults = {1000.0, static_cast< phred_t >( 'T' )};
    BOOST_CHECK_EQUAL_COLLECTIONS( expectedResults.begin(), expectedResults.end(), qualities.begin(), qualities.end() );
}

BOOST_AUTO_TEST_CASE( testShouldRetrieveCorrectQualitiesForInsertion )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "CATGA" );

    auto var1 = std::make_shared< Variant >( alignedSequence, Region( "1", 1, 1 ), "CAT" );
    // Check alignment correctly altered
    auto var2 = std::make_shared< Variant >( alignedSequence, Region( "1", 1, 2 ), "T" );

    const auto qualities = echidna::variant::getVariantsReadBaseQualities( 0, "STEFA", {var1, var2}, {} );

    const auto expectedResults = {static_cast< phred_t >( 'T' ), static_cast< phred_t >( 'A' )};
    BOOST_CHECK_EQUAL_COLLECTIONS( expectedResults.begin(), expectedResults.end(), qualities.begin(), qualities.end() );
}

BOOST_AUTO_TEST_CASE( testShouldRetrieveCorrectQualitiesWithBreakpointAtStart )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "CATGA" );

    auto bp = std::make_shared< Breakpoint >( "1", 0, true, "TAG" );
    // Check alignment correctly altered
    auto var = std::make_shared< Variant >( alignedSequence, Region( "1", 1, 2 ), "T" );

    const auto qualities = echidna::variant::getVariantsReadBaseQualities( 0, "STEFA", {var}, {bp} );

    const auto expectedResults = {static_cast< phred_t >( 'A' )};
    BOOST_CHECK_EQUAL_COLLECTIONS( expectedResults.begin(), expectedResults.end(), qualities.begin(), qualities.end() );
}

BOOST_AUTO_TEST_CASE( testShouldRetrieveCorrectQualitiesWithBreakpointNotAtStart )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 99, 105 ), "ACATGA" );

    auto bp = std::make_shared< Breakpoint >( "1", 105, false, "TAG" );
    // Check alignment correctly altered
    auto var = std::make_shared< Variant >( alignedSequence, Region( "1", 101, 102 ), "T" );

    const auto qualities = echidna::variant::getVariantsReadBaseQualities( 100, "STEFA", {var}, {bp} );

    const auto expectedResults = {static_cast< phred_t >( 'T' )};
    BOOST_CHECK_EQUAL_COLLECTIONS( expectedResults.begin(), expectedResults.end(), qualities.begin(), qualities.end() );
}

//-------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( testNormaliseVariantsOnStrandReturnsEmptyFromEmptyInput )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AGGGCACAGC" );
    const auto variants = echidna::variant::normaliseVariantsOnStrand( {}, *alignedSequence );

    BOOST_CHECK( variants.empty() );
}

BOOST_AUTO_TEST_CASE( testNormaliseVariantsCancelsInputs )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AGGGCACAGC" );

    auto del = std::make_shared< Variant >( alignedSequence, Region( "1", 1, 2 ), "" );
    auto ins = std::make_shared< Variant >( alignedSequence, Region( "1", 2, 2 ), "G" );

    const auto variants = echidna::variant::normaliseVariantsOnStrand( {del, ins}, *alignedSequence );

    BOOST_CHECK( variants.empty() );
}

BOOST_AUTO_TEST_CASE( testNormaliseVariantsOnStrandLeavesIsolatedSNPUnchanged )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AGGGCACAGC" );
    auto snp = std::make_shared< Variant >( alignedSequence, Region( "1", 2, 3 ), "A" );

    const auto variants = echidna::variant::normaliseVariantsOnStrand( {snp}, *alignedSequence );

    BOOST_CHECK_EQUAL( variants.size(), 1 );
    BOOST_CHECK( checkVariantInVector( std::vector< varPtr_t >( variants.begin(), variants.end() ),
                                       std::make_shared< Variant >( alignedSequence, Region( "1", 2, 3 ), "A" ) ) );
}

BOOST_AUTO_TEST_CASE( testNormaliseVariantsLeftAlignsInsertionToMinPos )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "GGGGGGGGGG" );
    auto var = std::make_shared< Variant >( alignedSequence, Region( "1", 3, 3 ), "G" );

    const auto variants = echidna::variant::normaliseVariantsOnStrand( {var}, *alignedSequence );

    BOOST_CHECK_EQUAL( variants.size(), 1 );
    BOOST_CHECK( checkVariantInVector( std::vector< varPtr_t >( variants.begin(), variants.end() ),
                                       std::make_shared< Variant >( alignedSequence, Region( "1", 0, 0 ), "G" ) ) );
}

BOOST_AUTO_TEST_CASE( testNormaliseVariantsLeftAlignsDeletionToMinPos )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "GGGGGGGGGG" );
    auto var = std::make_shared< Variant >( alignedSequence, Region( "1", 3, 4 ), "" );

    const auto variants = echidna::variant::normaliseVariantsOnStrand( {var}, *alignedSequence );

    BOOST_CHECK_EQUAL( variants.size(), 1 );
    BOOST_CHECK( checkVariantInVector( std::vector< varPtr_t >( variants.begin(), variants.end() ),
                                       std::make_shared< Variant >( alignedSequence, Region( "1", 0, 1 ), "" ) ) );
}

BOOST_AUTO_TEST_CASE( testNormaliseVariantsLeftAlignsDeletionToPositionOfSNP )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "GGGGGGGGGG" );

    auto snp = std::make_shared< Variant >( alignedSequence, Region( "1", 3, 4 ), "A" );
    auto var = std::make_shared< Variant >( alignedSequence, Region( "1", 5, 6 ), "" );

    const auto variants = echidna::variant::normaliseVariantsOnStrand( {snp, var}, *alignedSequence );

    BOOST_CHECK_EQUAL( variants.size(), 2 );
    BOOST_CHECK( checkVariantInVector( std::vector< varPtr_t >( variants.begin(), variants.end() ),
                                       std::make_shared< Variant >( alignedSequence, Region( "1", 0, 1 ), "" ) ) );
    BOOST_CHECK( checkVariantInVector( std::vector< varPtr_t >( variants.begin(), variants.end() ),
                                       std::make_shared< Variant >( alignedSequence, Region( "1", 4, 5 ), "A" ) ) );
}

BOOST_AUTO_TEST_CASE( testNormaliseVariantsLeftAlignsComplexCombinationOfInsertions )
{
    const auto alignedSequence =
        std::make_shared< ReferenceSequence >( Region( "1", 0, 29 ), "CATGATGATGATGATGATATATAAAAAAC" );

    auto ins1 = std::make_shared< Variant >( alignedSequence, Region( "1", 18, 18 ), "G" );
    auto ins2 = std::make_shared< Variant >( alignedSequence, Region( "1", 23, 23 ), "T" );
    auto ins3 = std::make_shared< Variant >( alignedSequence, Region( "1", 26, 26 ), "A" );

    const auto variants = echidna::variant::normaliseVariantsOnStrand( {ins1, ins2, ins3}, *alignedSequence );

    BOOST_CHECK_EQUAL( variants.size(), 1 );

    BOOST_CHECK( checkVariantInVector( std::vector< varPtr_t >( variants.begin(), variants.end() ),
                                       std::make_shared< Variant >( alignedSequence, Region( "1", 1, 1 ), "ATG" ) ) );
}

BOOST_AUTO_TEST_CASE( testNormalizationJoinsSNPsWithIndelsThenSplits )
{
    const auto alignedSequence =
        std::make_shared< ReferenceSequence >( Region( "1", 0, 29 ), "CATGATGATGATGATGATATATAAAAAAC" );

    auto ins = std::make_shared< Variant >( alignedSequence, Region( "1", 18, 18 ), "G" );
    auto snp = std::make_shared< Variant >( alignedSequence, Region( "1", 18, 19 ), "G" );

    const auto variants = echidna::variant::normaliseVariantsOnStrand( {ins, snp}, *alignedSequence );

    BOOST_CHECK_EQUAL( variants.size(), 2 );

    BOOST_CHECK( checkVariantInVector( std::vector< varPtr_t >( variants.begin(), variants.end() ),
                                       std::make_shared< Variant >( alignedSequence, Region( "1", 18, 18 ), "G" ) ) );
    BOOST_CHECK( checkVariantInVector( std::vector< varPtr_t >( variants.begin(), variants.end() ),
                                       std::make_shared< Variant >( alignedSequence, Region( "1", 18, 19 ), "G" ) ) );
}

BOOST_AUTO_TEST_CASE( testNormaliseVariantsLeftAlignsInsertionToPositionOfSNPJoinsAndTrimms )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "GGGGGGGGGG" );

    auto snp = std::make_shared< Variant >( alignedSequence, Region( "1", 3, 4 ), "A" );
    auto var = std::make_shared< Variant >( alignedSequence, Region( "1", 8, 8 ), "G" );

    const auto variants = echidna::variant::normaliseVariantsOnStrand( {snp, var}, *alignedSequence );

    BOOST_CHECK_EQUAL( variants.size(), 1 );
    BOOST_CHECK( checkVariantInVector( std::vector< varPtr_t >( variants.begin(), variants.end() ),
                                       std::make_shared< Variant >( alignedSequence, Region( "1", 3, 3 ), "A" ) ) );
}

BOOST_AUTO_TEST_CASE( testNormaliseVariantsLeftAlignsDeletionJoinsInsertionsTogether )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "TGGGGGGGGG" );

    auto var1 = std::make_shared< Variant >( alignedSequence, Region( "1", 3, 4 ), "" );
    auto var2 = std::make_shared< Variant >( alignedSequence, Region( "1", 5, 6 ), "" );

    const auto variants = echidna::variant::normaliseVariantsOnStrand( {var1, var2}, *alignedSequence );

    BOOST_CHECK_EQUAL( variants.size(), 1 );
    BOOST_CHECK( checkVariantInVector( std::vector< varPtr_t >( variants.begin(), variants.end() ),
                                       std::make_shared< Variant >( alignedSequence, Region( "1", 1, 3 ), "" ) ) );
}

BOOST_AUTO_TEST_CASE( testNormaliseVariantsLeftAlignsInsertionJoinsInsertionsTogether )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "TGGGGGGGGG" );

    auto var1 = std::make_shared< Variant >( alignedSequence, Region( "1", 6, 6 ), "G" );
    auto var2 = std::make_shared< Variant >( alignedSequence, Region( "1", 8, 8 ), "G" );

    const auto variants = echidna::variant::normaliseVariantsOnStrand( {var1, var2}, *alignedSequence );

    BOOST_CHECK_EQUAL( variants.size(), 1 );
    BOOST_CHECK( checkVariantInVector( std::vector< varPtr_t >( variants.begin(), variants.end() ),
                                       std::make_shared< Variant >( alignedSequence, Region( "1", 1, 1 ), "GG" ) ) );
}

BOOST_AUTO_TEST_CASE( testNormalisationCancelsDeletionAndIndel )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "TGGGGGGGGG" );
    auto var1 = std::make_shared< Variant >( alignedSequence, Region( "1", 6, 8 ), "" );
    auto var2 = std::make_shared< Variant >( alignedSequence, Region( "1", 8, 8 ), "GG" );

    const auto variants = echidna::variant::normaliseVariantsOnStrand( {var1, var2}, *alignedSequence );

    BOOST_CHECK_EQUAL( variants.size(), 0 );
}

BOOST_AUTO_TEST_CASE( testNormalisationCancelsDeletionAndIndelToProduceSNPs )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "ATCGTCTAGT" );
    auto var1 = std::make_shared< Variant >( alignedSequence, Region( "1", 6, 8 ), "" );
    auto var2 = std::make_shared< Variant >( alignedSequence, Region( "1", 8, 8 ), "GG" );

    const auto variants = echidna::variant::normaliseVariantsOnStrand( {var1, var2}, *alignedSequence );

    BOOST_CHECK_EQUAL( variants.size(), 2 );
    BOOST_CHECK( checkVariantInVector( std::vector< varPtr_t >( variants.begin(), variants.end() ),
                                       std::make_shared< Variant >( alignedSequence, Region( "1", 6, 7 ), "G" ) ) );
    BOOST_CHECK( checkVariantInVector( std::vector< varPtr_t >( variants.begin(), variants.end() ),
                                       std::make_shared< Variant >( alignedSequence, Region( "1", 7, 8 ), "G" ) ) );
}

BOOST_AUTO_TEST_CASE( testVariantGeneratorSingleReadWithNoIndels )
{
    int64_t startPos = 500L;
    const auto alignedSequence = std::make_shared< ReferenceSequence >(
        Region( "1", startPos, startPos + 100 ),
        "AGGGCACAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACGGGAAGGA" );

    auto mutSeqStr = alignedSequence->sequence().str();

    mutSeqStr[0] = 'C';
    mutSeqStr[1] = 'T';
    mutSeqStr[11] = 'C';
    mutSeqStr[22] = 'T';
    mutSeqStr[52] = 'T';

    auto read = std::make_shared< Read >( mutSeqStr, std::string( mutSeqStr.size(), 'Q' ), "-", Cigar( "100M" ), 0,
                                          startPos, 0, 100, 0, 200, 200, alignedSequence );

    ReadDataset readDataset( {"sample"}, alignedSequence->region() );
    readDataset.insertRead( {"sample"}, read );

    VariantGenerator variantGenerator( alignedSequence, 10, 0 );

    VariantContainer variantContainer = variantGenerator.generateVariantsFromReads( readDataset.getAllReads( 0 ) );

    const auto variants = variantContainer.getVariants();
    const std::vector< varPtr_t > vecVariants( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( vecVariants.size(), 5 );

    BOOST_CHECK_EQUAL( *vecVariants[0], Variant( alignedSequence, Region( "1", startPos + 0, startPos + 1 ), "C" ) );
    BOOST_CHECK_EQUAL( *vecVariants[1], Variant( alignedSequence, Region( "1", startPos + 1, startPos + 2 ), "T" ) );
    BOOST_CHECK_EQUAL( *vecVariants[2], Variant( alignedSequence, Region( "1", startPos + 11, startPos + 12 ), "C" ) );
    BOOST_CHECK_EQUAL( *vecVariants[3], Variant( alignedSequence, Region( "1", startPos + 22, startPos + 23 ), "T" ) );
    BOOST_CHECK_EQUAL( *vecVariants[4], Variant( alignedSequence, Region( "1", startPos + 52, startPos + 53 ), "T" ) );
}

BOOST_AUTO_TEST_CASE( testShouldNotCountSupportOfLowQualityVariant )
{
    int64_t startPos = 500L;
    const auto alignedSequence = std::make_shared< ReferenceSequence >(
        Region( "1", startPos, startPos + 100 ),
        "AGGGCACAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACGGGAAGGA" );

    const char minBaseQuality = 10;
    auto mutSeqStr = alignedSequence->sequence().str();
    auto qualString = std::string( mutSeqStr.size(), 'Q' );

    mutSeqStr[0] = 'C';
    qualString[0] = minBaseQuality;

    mutSeqStr[1] = 'T';
    qualString[1] = minBaseQuality - 1;

    auto read = std::make_shared< Read >( mutSeqStr, qualString, "-", Cigar( "100M" ), 0, startPos, 0, 100, 0, 200, 200,
                                          alignedSequence );

    ReadDataset readDataset( {"sample"}, alignedSequence->region() );
    readDataset.insertRead( {"sample"}, read );

    VariantGenerator variantGenerator( alignedSequence, minBaseQuality, 0 );

    VariantContainer variantContainer = variantGenerator.generateVariantsFromReads( readDataset.getAllReads( 0 ) );

    const auto variants = variantContainer.getVariants();
    const std::vector< varPtr_t > vecVariants( variants.begin(), variants.end() );

    BOOST_REQUIRE_EQUAL( vecVariants.size(), 2 );
    BOOST_CHECK_EQUAL( *vecVariants[0], Variant( alignedSequence, Region( "1", startPos + 0, startPos + 1 ), "C" ) );
    BOOST_CHECK_EQUAL( *vecVariants[1], Variant( alignedSequence, Region( "1", startPos + 1, startPos + 2 ), "T" ) );

    BOOST_CHECK_EQUAL( variantContainer.totalReadsSupportingVariant( vecVariants[0] ), 1 );
    BOOST_CHECK_EQUAL( variantContainer.totalReadsSupportingVariant( vecVariants[1] ), 0 );
}

BOOST_AUTO_TEST_CASE( testShouldNotReadWithLowSupport )
{
    int64_t startPos = 500L;
    const auto alignedSequence = std::make_shared< ReferenceSequence >(
        Region( "1", startPos, startPos + 100 ),
        "AGGGCACAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACGGGAAGGA" );

    const uint8_t minMappingQuality = 10;
    auto mutSeqStr = alignedSequence->sequence().str();
    auto qualString = std::string( mutSeqStr.size(), 'Q' );

    mutSeqStr[0] = 'C';

    auto read1 = std::make_shared< Read >( mutSeqStr, qualString, "-", Cigar( "100M" ), 0, startPos, 0,
                                           minMappingQuality - 1, 0, 200, 200, alignedSequence );

    auto read2 = std::make_shared< Read >( mutSeqStr, qualString, "-", Cigar( "100M" ), 0, startPos, 0,
                                           minMappingQuality, 0, 200, 200, alignedSequence );

    ReadDataset readDataset( {"sample"}, alignedSequence->region() );
    readDataset.insertRead( {"sample"}, read1 );
    readDataset.insertRead( {"sample"}, read2 );

    VariantGenerator variantGenerator( alignedSequence, 0, minMappingQuality );

    VariantContainer variantContainer = variantGenerator.generateVariantsFromReads( readDataset.getAllReads( 0 ) );

    const auto variants = variantContainer.getVariants();
    const std::vector< varPtr_t > vecVariants( variants.begin(), variants.end() );

    BOOST_REQUIRE_EQUAL( vecVariants.size(), 1 );
    BOOST_CHECK_EQUAL( *vecVariants[0], Variant( alignedSequence, Region( "1", startPos + 0, startPos + 1 ), "C" ) );

    BOOST_CHECK_EQUAL( variantContainer.totalReadsSupportingVariant( vecVariants[0] ), 1 );
}

BOOST_AUTO_TEST_CASE( testVariantGeneratorMulipleReadsWithNoIndels )
{
    int64_t startPos = 500L;
    const auto alignedSequence = std::make_shared< ReferenceSequence >(
        Region( "1", startPos, startPos + 100 ),
        "AGGGCACAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACGGGAAGGA" );

    auto mutSeqStr1 = alignedSequence->sequence().str();
    auto mutSeqStr2 = alignedSequence->sequence().str();

    mutSeqStr1[0] = 'C';
    mutSeqStr2[1] = 'T';

    mutSeqStr1[11] = 'C';
    mutSeqStr2[11] = 'C';

    mutSeqStr1[22] = 'T';
    mutSeqStr2[52] = 'T';

    auto read1 = std::make_shared< Read >( mutSeqStr1, std::string( mutSeqStr1.size(), 'Q' ), "-", Cigar( "100M" ), 0,
                                           startPos, 0, 100, 0, 200, 200, alignedSequence );

    auto read2 = std::make_shared< Read >( mutSeqStr2, std::string( mutSeqStr2.size(), 'R' ), "-", Cigar( "100M" ), 0,
                                           startPos, 0, 100, 0, 200, 200, alignedSequence );

    ReadDataset readDataset( {"sample"}, alignedSequence->region() );
    readDataset.insertRead( {"sample"}, read1 );
    readDataset.insertRead( {"sample"}, read2 );

    VariantGenerator variantGenerator( alignedSequence, 10, 0 );

    VariantContainer variantContainer = variantGenerator.generateVariantsFromReads( readDataset.getAllReads( 0 ) );

    const auto variants = variantContainer.getVariants();
    const std::vector< varPtr_t > vecVariants( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( vecVariants.size(), 5 );

    BOOST_CHECK_EQUAL( *vecVariants[0], Variant( alignedSequence, Region( "1", startPos + 0, startPos + 1 ), "C" ) );
    BOOST_CHECK_EQUAL( *vecVariants[1], Variant( alignedSequence, Region( "1", startPos + 1, startPos + 2 ), "T" ) );
    BOOST_CHECK_EQUAL( *vecVariants[2], Variant( alignedSequence, Region( "1", startPos + 11, startPos + 12 ), "C" ) );
    BOOST_CHECK_EQUAL( *vecVariants[3], Variant( alignedSequence, Region( "1", startPos + 22, startPos + 23 ), "T" ) );
    BOOST_CHECK_EQUAL( *vecVariants[4], Variant( alignedSequence, Region( "1", startPos + 52, startPos + 53 ), "T" ) );
}

BOOST_AUTO_TEST_CASE( testVariantGeneratorSingleReadWithDeletionStartingBeforeBlockRegion )
{
    const auto alignedSequence = std::make_shared< ReferenceSequence >(
        echidna::caller::Region( "1", 499, 600 ),
        "TAGGGCACAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACGGGAAGGA" );

    int64_t startPos = 500L;
    auto mutSeqStr =
        alignedSequence->subseq( Region( "1", startPos, startPos + 100 ) ).sequence().substr( 2, 98 ).str();
    mutSeqStr[11 - 2] = 'C';
    mutSeqStr[22 - 2] = 'T';
    mutSeqStr[52 - 2] = 'T';  // Should not find this snp as beyond the region boundary.
    echidna::utils::BasePairSequence mutSeq = mutSeqStr;

    auto read = std::make_shared< Read >( mutSeqStr, std::string( mutSeqStr.size(), 'Q' ), "-", Cigar( "3D98M" ), 0,
                                          startPos - 1, 0, 100, 200, 200, 0, alignedSequence );

    ReadDataset readDataset( {"sample"}, Region( "1", startPos, startPos ) );
    readDataset.insertRead( {"sample"}, read );

    VariantGenerator variantGenerator( alignedSequence, 10, 0 );

    VariantContainer variantContainer = variantGenerator.generateVariantsFromReads( readDataset.getAllReads( 0 ) );

    auto variants = variantContainer.getVariants();
    std::vector< varPtr_t > vecVariants( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( vecVariants.size(), 3 );

    BOOST_CHECK( checkVariantInVector(
        vecVariants,
        std::make_shared< Variant >( alignedSequence, Region( "1", startPos + 11, startPos + 12 ), "C" ) ) );
    BOOST_CHECK( checkVariantInVector(
        vecVariants,
        std::make_shared< Variant >( alignedSequence, Region( "1", startPos + 22, startPos + 23 ), "T" ) ) );

    BOOST_CHECK_EQUAL( *vecVariants[2], Variant( alignedSequence, Region( "1", startPos + 52, startPos + 53 ), "T" ) );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testVariantGeneratorJoinsSnpWithIndel )
{

    const auto alignedSequence = std::make_shared< ReferenceSequence >(
        Region( "1", 0L, 100L ),
        "NGGGCACAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAGACCCCCTGCCAGGGAGACCCCAGGCGCCTGAATGGCCACGGGAAGGA" );

    auto mutSeqWithOutIndelsStr = alignedSequence->sequence().str();
    mutSeqWithOutIndelsStr[29] = 'T';  // Acts as a block to left-alignment of cigarFlags::DELETION on right.
    echidna::utils::BasePairSequence mutSeqWithOutIndels = mutSeqWithOutIndelsStr;

    auto mutSeq = mutSeqWithOutIndels.substr( 7, 25 ) + mutSeqWithOutIndels.substr( 7 + 25 + 1, 65 );

    const auto read = std::make_shared< Read >( mutSeq, std::string( mutSeq.size(), 'Q' ), "test_read",
                                                Cigar( "7D25M1D65M2D" ), 0, 0, 0, 100, 200, 200, 0, alignedSequence );

    ReadDataset readDataset( {"sample"}, Region( "1", 0, 100 ) );
    readDataset.insertRead( {"sample"}, read );

    // Get padded refSeq.
    auto variants = VariantGenerator( alignedSequence, 10, 0 )
                        .generateVariantsFromReads( readDataset.getAllReads( 0 ) )
                        .getVariants();

    std::vector< varPtr_t > vecVariants( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( vecVariants.size(), 2 );

    BOOST_CHECK( checkVariantInVector( vecVariants,
                                       std::make_shared< Variant >( alignedSequence, Region( "1", 29, 30 ), "" ) ) );

    BOOST_CHECK( checkVariantInVector( vecVariants,
                                       std::make_shared< Variant >( alignedSequence, Region( "1", 30, 31 ), "T" ) ) );
    for ( const auto & var : vecVariants )
    {
        const auto readsFromVar = var->getReads();
        std::vector< readPtr_t > expected = {read};
        BOOST_CHECK_EQUAL_COLLECTIONS( expected.begin(), expected.end(), readsFromVar.begin(), readsFromVar.end() );
    }
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( shouldNotCallIndelsAtReadStartEnd )
{
    std::string varGenRef( 7, 'A' );
    const auto alignedSequence =
        std::make_shared< ReferenceSequence >( echidna::caller::Region( "1", -1L, 7L ), "A" + varGenRef );

    auto read = std::make_shared< Read >( "T" + std::string( 6, 'A' ) + "TT", std::string( 9, 'Q' ), "",
                                          Cigar( "1I6M2I" ), 0, 0, 0, 0, 100, 200, 200, alignedSequence );

    ReadDataset readDataset( {"sample"}, Region( "1", 0, 100 ) );
    readDataset.insertRead( {"sample"}, read );

    VariantGenerator variantGenerator( alignedSequence, 10, 0 );
    auto variantContainer = variantGenerator.generateVariantsFromReads( readDataset.getAllReads( 0 ) );

    auto variants = variantContainer.getVariants();

    BOOST_CHECK( variants.empty() );
}
