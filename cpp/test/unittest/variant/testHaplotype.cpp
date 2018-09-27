// All content Copyright (C) 2018 Genomics plc
#include "variant/haplotype.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <string>

#include "variant/type/variant.hpp"

using namespace echidna::variant;
using echidna::variant::variantSet_t;
using echidna::caller::Region;
using echidna::caller::SetRegions;
using echidna::utils::ReferenceSequence;

BOOST_AUTO_TEST_CASE( testBuildingHaplotypeSequenceWithSubDeletion )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "TATCG" );
    const auto var1 = std::make_shared< Variant >( referenceSequence, Region( "1", 2, 3 ), "" );

    BOOST_CHECK_EQUAL( Haplotype::buildHaplotypeSequence( referenceSequence, referenceSequence->region(), {var1} ),
                       "TACG" );
}

BOOST_AUTO_TEST_CASE( testBuildingHaplotypeSequenceFullDeletion )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "TAAAA" );
    const auto var1 = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "" );

    BOOST_CHECK_EQUAL( Haplotype::buildHaplotypeSequence( referenceSequence, referenceSequence->region(), {var1} ),
                       "" );
}

BOOST_AUTO_TEST_CASE( testBuildingHaplotypeSequenceFullMNP )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "TAAAA" );
    const auto var1 = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "ATTTT" );

    BOOST_CHECK_EQUAL( Haplotype::buildHaplotypeSequence( referenceSequence, referenceSequence->region(), {var1} ),
                       "ATTTT" );
}

BOOST_AUTO_TEST_CASE( testBuildingHaplotypeSequenceInsertionAtStart )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "TAAAA" );
    const auto var1 = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 0 ), "AGGGG" );

    BOOST_CHECK_EQUAL( Haplotype::buildHaplotypeSequence( referenceSequence, referenceSequence->region(), {var1} ),
                       "AGGGGTAAAA" );
}

BOOST_AUTO_TEST_CASE( testBuildingHaplotypeSequenceInsertionAtEnd )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "TAAAA" );
    const auto var1 = std::make_shared< Variant >( referenceSequence, Region( "1", 5, 5 ), "AGGGG" );

    BOOST_CHECK_EQUAL( Haplotype::buildHaplotypeSequence( referenceSequence, referenceSequence->region(), {var1} ),
                       "TAAAAAGGGG" );
}

BOOST_AUTO_TEST_CASE( testIsValidVariantCombinationShouldCheckForOverlap )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "TAAAA" );
    auto var1 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "T" );
    auto var2 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "G" );
    const variantSet_t variants = {var1, var2};

    BOOST_CHECK( not isValidVariantCombination( variants.cbegin(), variants.cend(), referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( testIsValidVariantCombinationShouldReturnFalseForTwoInsertionsAtSameLocation )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "TAAAA" );
    auto var1 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 1 ), "AT" );
    auto var2 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 1 ), "AG" );
    const variantSet_t variants = {var1, var2};

    BOOST_CHECK( not isValidVariantCombination( variants.cbegin(), variants.cend(), referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( testIsValidVariantCombinationShouldCheckForDifferentRepresentationsOfIndels )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "TAAAA" );
    auto var1 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "" );
    auto var2 = std::make_shared< Variant >( referenceSequence, Region( "1", 3, 4 ), "" );
    const variantSet_t variants = {var1, var2};

    BOOST_CHECK( not isValidVariantCombination( variants.cbegin(), variants.cend(), referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( testIsValidVariantCombinationShouldNotCheckIndelsSeparatedByOtherVariants )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "TAAAA" );
    auto var1 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "" );
    auto var2 = std::make_shared< Variant >( referenceSequence, Region( "1", 2, 3 ), "G" );
    auto var3 = std::make_shared< Variant >( referenceSequence, Region( "1", 3, 4 ), "" );
    const variantSet_t variants = {var1, var2, var3};

    BOOST_CHECK( isValidVariantCombination( variants.cbegin(), variants.cend(), referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( testIsValidVariantCombinationShouldFailWhenReferenceDoesntCoverFirstVariant )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "TAAAA" );
    auto var1 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "" );
    auto var2 = std::make_shared< Variant >( referenceSequence, Region( "1", 3, 4 ), "" );
    const variantSet_t variants = {var1, var2};

    auto incompleteReference =
        std::make_shared< ReferenceSequence >( referenceSequence->subseq( Region( "1", 2, 5 ) ) );

    BOOST_CHECK_THROW( isValidVariantCombination( variants.cbegin(), variants.cend(), incompleteReference ),
                       echidna::utils::echidna_exception );
}

BOOST_AUTO_TEST_CASE( testIsValidVariantCombinationShouldFailWhenReferenceDoesntCoverSecondVariant )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 5 ), "TAAAA" );
    auto var1 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "" );
    auto var2 = std::make_shared< Variant >( referenceSequence, Region( "1", 3, 4 ), "" );
    const variantSet_t variants = {var1, var2};

    auto incompleteReference =
        std::make_shared< ReferenceSequence >( referenceSequence->subseq( Region( "1", 0, 2 ) ) );

    BOOST_CHECK_THROW( isValidVariantCombination( variants.cbegin(), variants.cend(), incompleteReference ),
                       echidna::utils::echidna_exception );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testSingleSNPHap )
{
    std::string chrom = "1";

    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( chrom, 1000007, 1000013 ), "GCCTCA" );
    const varPtr_t var = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000010, 1000011 ), "A" );
    const variantSet_t variants = {var};

    const Haplotype theHap( referenceSequence, referenceSequence->region(), variants, 0 );

    const std::string mutSeq( "GCCACA" );

    BOOST_CHECK_EQUAL( theHap.paddedSequences().front(), mutSeq );

    BOOST_CHECK( theHap.isReference( Region( chrom, referenceSequence->start(), var->start() ) ) );
    BOOST_CHECK( not theHap.isReference( var->region() ) );
    BOOST_CHECK( theHap.isReference( Region( chrom, var->end(), referenceSequence->end() ) ) );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testSingleInsertionHap )
{
    using namespace echidna::variant;

    std::string chrom = "1";

    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( chrom, 1000007, 1000013 ), "GCCTCA" );
    const varPtr_t var =
        std::make_shared< Variant >( referenceSequence, Region( chrom, 1000011, 1000011 ), "ACTGACTG" );
    const variantSet_t variants = {var};

    const Haplotype theHap( referenceSequence, referenceSequence->region(), variants, 0 );
    const std::string mutSeq( "GCCTACTGACTGCA" );

    BOOST_CHECK_EQUAL( theHap.paddedSequences().front(), mutSeq );

    BOOST_CHECK( theHap.isReference( Region( chrom, referenceSequence->start(), var->start() ) ) );
    BOOST_CHECK( not theHap.isReference( var->region() ) );
    BOOST_CHECK( theHap.isReference( Region( chrom, var->end(), referenceSequence->end() ) ) );
}

////-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testSingleDeletionHap )
{
    using namespace echidna::variant;

    std::string chrom = "1";

    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( chrom, 1000007, 1000017 ), "GCCTCCACCA" );
    const varPtr_t var = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000011, 1000014 ), "" );
    const variantSet_t variants = {var};
    const Haplotype theHap( referenceSequence, referenceSequence->region(), variants, 0 );

    const std::string mutSeq( "GCCTCCA" );

    BOOST_CHECK_EQUAL( theHap.paddedSequences().front(), mutSeq );

    BOOST_CHECK( theHap.isReference( Region( chrom, referenceSequence->start(), var->start() ) ) );
    BOOST_CHECK( not theHap.isReference( var->region() ) );
    BOOST_CHECK( theHap.isReference( Region( chrom, var->end(), referenceSequence->end() ) ) );
}

////-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testSingleMNPHap )
{
    using namespace echidna::variant;

    std::string chrom = "1";

    auto reference = std::make_shared< ReferenceSequence >( Region( chrom, 1000007, 1000014 ), "GCCTCAC" );
    const varPtr_t var = std::make_shared< Variant >( reference, Region( chrom, 1000010, 1000013 ), "ACT" );
    const variantSet_t variants = {var};
    const Haplotype theHap( reference, reference->region(), variants, 0 );

    const std::string mutSeq( "GCCACTC" );

    BOOST_CHECK_EQUAL( theHap.paddedSequences().front(), mutSeq );

    BOOST_CHECK( theHap.containsVariant( var ) );

    BOOST_CHECK( theHap.isReference( Region( chrom, reference->start(), var->start() ) ) );
    BOOST_CHECK( not theHap.isReference( var->region() ) );
    BOOST_CHECK( not theHap.isReference( Region( chrom, 1000010, 1000011 ) ) );
    BOOST_CHECK( theHap.isReference( Region( chrom, 1000011, 1000012 ) ) );
    BOOST_CHECK( not theHap.isReference( Region( chrom, 1000012, 1000013 ) ) );
    BOOST_CHECK( theHap.isReference( Region( chrom, var->end(), reference->end() ) ) );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testMultiSNPHap )
{
    using namespace echidna::variant;

    std::string chrom = "1";

    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( chrom, 1000007, 1000018 ), "GCCTCACCCAG" );
    const varPtr_t var1 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000010, 1000011 ), "A" );
    const varPtr_t var2 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000012, 1000013 ), "G" );
    const varPtr_t var3 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000015, 1000016 ), "T" );

    const variantSet_t variants = {var1, var2, var3};

    const Haplotype theHap( referenceSequence, referenceSequence->region(), variants, 0 );
    const std::string mutSeq( "GCCACGCCTAG" );

    BOOST_CHECK_EQUAL( theHap.paddedSequences().front(), mutSeq );

    BOOST_CHECK( theHap.containsVariant( var1 ) );
    BOOST_CHECK( theHap.containsVariant( var2 ) );
    BOOST_CHECK( theHap.containsVariant( var3 ) );

    BOOST_CHECK( theHap.isReference( Region( chrom, referenceSequence->start(), var1->start() ) ) );

    BOOST_CHECK( not theHap.isReference( var1->region() ) );
    BOOST_CHECK( theHap.isReference( Region( chrom, 1000011, 1000012 ) ) );
    BOOST_CHECK( not theHap.isReference( var2->region() ) );
    BOOST_CHECK( not theHap.isReference( var3->region() ) );

    BOOST_CHECK( theHap.isReference( Region( chrom, var3->end(), referenceSequence->end() ) ) );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testMultiInsertionHap )
{
    using namespace echidna::variant;

    std::string chrom = "1";

    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( chrom, 1000007, 1000018 ), "GCCTCACCCAG" );
    const varPtr_t var1 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000011, 1000011 ), "A" );
    const varPtr_t var2 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000013, 1000013 ), "A" );
    const varPtr_t var3 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000016, 1000016 ), "A" );

    const variantSet_t variants = {var1, var2, var3};

    const Haplotype theHap( referenceSequence, referenceSequence->region(), variants, 0 );

    const std::string mutSeq( "GCCTACAACCCAAG" );

    BOOST_CHECK_EQUAL( theHap.paddedSequences().front(), mutSeq );

    BOOST_CHECK( theHap.containsVariant( var1 ) );
    BOOST_CHECK( theHap.containsVariant( var2 ) );
    BOOST_CHECK( theHap.containsVariant( var3 ) );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testMultiDeletionHap )
{
    using namespace echidna::variant;

    std::string chrom = "1";

    auto referenceSequence =
        std::make_shared< ReferenceSequence >( Region( chrom, 1000007, 1000024 ), "GCCTCACCCAGGAAAGC" );
    const varPtr_t var1 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000011, 1000012 ), "" );
    const varPtr_t var2 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000014, 1000016 ), "" );
    const varPtr_t var3 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000019, 1000022 ), "" );

    const variantSet_t variants = {var1, var2, var3};

    const Haplotype theHap( referenceSequence, referenceSequence->region(), variants, 0 );

    const std::string mutSeq( "GCCTACAGGGC" );

    BOOST_CHECK_EQUAL( theHap.paddedSequences().front(), mutSeq );

    BOOST_CHECK( theHap.containsVariant( var1 ) );
    BOOST_CHECK( theHap.containsVariant( var2 ) );
    BOOST_CHECK( theHap.containsVariant( var3 ) );
}

////-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testMNPShouldBeContainedInHaplotype )
{
    std::string chrom = "1";

    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( chrom, 1000007, 1000018 ), "GCCTCACCCAG" );

    const varPtr_t var1 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000010, 1000011 ), "A" );
    const varPtr_t var2 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000012, 1000013 ), "G" );
    const varPtr_t var3 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000015, 1000016 ), "T" );
    const varPtr_t mnp = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000010, 1000016 ), "ACGCCT" );

    const variantSet_t variants = {var1, var2, var3};

    Haplotype theHap( referenceSequence, referenceSequence->region(), variants, 0 );
    const std::string mutSeq( "GCCACGCCTAG" );

    BOOST_CHECK_EQUAL( theHap.paddedSequences().front(), mutSeq );

    BOOST_CHECK( theHap.containsVariant( var1 ) );
    BOOST_CHECK( theHap.containsVariant( var2 ) );
    BOOST_CHECK( theHap.containsVariant( var3 ) );

    auto mnps = theHap.withMNPs();
    BOOST_CHECK_EQUAL( mnps.size(), 1 );

    BOOST_CHECK( not theHap.containsVariant( var1 ) );
    BOOST_CHECK( not theHap.containsVariant( var2 ) );
    BOOST_CHECK( not theHap.containsVariant( var3 ) );
    BOOST_CHECK( theHap.containsVariant( mnp ) );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testAdjacentSnpInsertion )
{
    using namespace echidna::variant;

    std::string chrom = "1";

    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( chrom, 1000007, 1000013 ), "GCCTCA" );
    const varPtr_t var1 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000010, 1000011 ), "A" );
    const varPtr_t var2 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000011, 1000011 ), "CCC" );

    const variantSet_t variants = {var1, var2};

    const Haplotype theHap( referenceSequence, referenceSequence->region(), variants, 0 );

    const std::string mutSeq( "GCCACCCCA" );

    BOOST_CHECK_EQUAL( theHap.paddedSequences().front(), mutSeq );
    BOOST_CHECK( theHap.containsVariant( var1 ) );
    BOOST_CHECK( theHap.containsVariant( var2 ) );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testAdjacentSnpDeletion )
{
    using namespace echidna::variant;

    std::string chrom = "1";

    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( chrom, 1000007, 1000016 ), "GCCTCACCC" );
    const varPtr_t var1 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000010, 1000011 ), "A" );
    const varPtr_t var2 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000012, 1000014 ), "" );

    const variantSet_t variants = {var1, var2};

    const Haplotype theHap( referenceSequence, referenceSequence->region(), variants, 0 );

    const std::string mutSeq( "GCCACCC" );

    BOOST_CHECK_EQUAL( theHap.paddedSequences().front(), mutSeq );
    BOOST_CHECK( theHap.containsVariant( var1 ) );
    BOOST_CHECK( theHap.containsVariant( var2 ) );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testInsertionsAtSameLocationShouldRaiseAsInvalidCombination )
{
    using namespace echidna::variant;

    std::string chrom = "20";

    auto referenceSequence =
        std::make_shared< ReferenceSequence >( Region( chrom, 106701, 106721 ), "CCTTTTTTTTTTTTTTTTGA" );
    const varPtr_t var1 = std::make_shared< Variant >( referenceSequence, Region( chrom, 106703, 106703 ), "T" );
    const varPtr_t var2 = std::make_shared< Variant >( referenceSequence, Region( chrom, 106703, 106703 ), "TTT" );

    const variantSet_t variants = {var1, var2};

    BOOST_CHECK_THROW( Haplotype( referenceSequence, referenceSequence->region(), variants, 0 ), std::runtime_error );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testInsertionAndDeletionAtSameLocationInHomopolymerRegion )
{
    using namespace echidna::variant;

    std::string chrom = "20";

    auto referenceSequence =
        std::make_shared< ReferenceSequence >( Region( chrom, 106701, 106721 ), "CCTTTTTTTTTTTTTTTTGA" );
    const varPtr_t var1 = std::make_shared< Variant >( referenceSequence, Region( chrom, 106703, 106703 ), "TTT" );
    const varPtr_t var2 = std::make_shared< Variant >( referenceSequence, Region( chrom, 106703, 106704 ), "" );

    varPtrComp comp;

    const variantSet_t variants = {var1, var2};
    const Haplotype theHap( referenceSequence, referenceSequence->region(), variants, 0 );

    const std::string mutSeq( "CCTTTTTTTTTTTTTTTTTTGA" );

    BOOST_CHECK( comp( var1, var2 ) );
    BOOST_CHECK_EQUAL( theHap.paddedSequences().front(), mutSeq );
    BOOST_CHECK( theHap.containsVariant( var1 ) );
    BOOST_CHECK( theHap.containsVariant( var2 ) );
}

BOOST_AUTO_TEST_CASE( testHaplotypeConstructionWithSubRegions )
{
    using namespace echidna::variant;

    std::string chrom = "20";

    auto referenceSequence =
        std::make_shared< ReferenceSequence >( Region( chrom, 106701, 106722 ), "CCTTTTTTTTTTTTTTTTGAT" );
    const auto var1 = std::make_shared< Variant >( referenceSequence, Region( chrom, 106703, 106703 ), "TTT" );
    const auto var2 = std::make_shared< Variant >( referenceSequence, Region( chrom, 106720, 106721 ), "" );

    const variantSet_t variants = {var1, var2};
    const auto padding = 0;

    SetRegions regions;
    regions.insert( Region( chrom, 106701, 106703 ) );
    regions.insert( Region( chrom, 106711, 106713 ) );
    regions.insert( Region( chrom, 106719, 106721 ) );

    const Haplotype theHap( referenceSequence, regions, variants, padding );

    const auto & sequences = theHap.paddedSequences();
    BOOST_REQUIRE_EQUAL( sequences.size(), 3 );
    BOOST_CHECK_EQUAL( "CCTTT", sequences[0] );
    BOOST_CHECK_EQUAL( "TT", sequences[1] );
    BOOST_CHECK_EQUAL( "G", sequences[2] );

    const auto & hapRegions = theHap.regions();
    BOOST_REQUIRE_EQUAL( hapRegions, regions );
}

BOOST_AUTO_TEST_CASE( testShouldJoinRegionsWhichHaveLargeVariantSpanning )
{
    using namespace echidna::variant;

    std::string chrom = "20";

    auto referenceSequence =
        std::make_shared< ReferenceSequence >( Region( chrom, 106701, 106722 ), "CCTTTTTTTTTTTTTTTTGAT" );
    const auto var1 = std::make_shared< Variant >( referenceSequence, Region( chrom, 106703, 106711 ), "TTT" );
    const auto var2 = std::make_shared< Variant >( referenceSequence, Region( chrom, 106720, 106721 ), "" );

    const variantSet_t variants = {var1, var2};
    const auto padding = 0;

    SetRegions regions;
    regions.insert( Region( chrom, 106701, 106703 ) );
    regions.insert( Region( chrom, 106711, 106713 ) );
    regions.insert( Region( chrom, 106719, 106721 ) );

    const Haplotype theHap( referenceSequence, regions, variants, padding );

    const auto & sequences = theHap.paddedSequences();
    BOOST_REQUIRE_EQUAL( sequences.size(), 2 );
    BOOST_CHECK_EQUAL( "CCTTTTT", sequences[0] );
    BOOST_CHECK_EQUAL( "G", sequences[1] );

    const auto & hapRegions = theHap.regions();
    SetRegions expectedRegions;
    expectedRegions.insert( Region( chrom, 106701, 106713 ) );
    expectedRegions.insert( Region( chrom, 106719, 106721 ) );
    BOOST_REQUIRE_EQUAL( hapRegions, expectedRegions );
}

BOOST_AUTO_TEST_CASE( testHaplotypeConstructionShouldFillGapsBetweenSubRegionsAndPadTheSequences )
{
    using namespace echidna::variant;

    std::string chrom = "20";

    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( chrom, 106690, 106730 ),
                                                                    "AAAAAAAAAAACCTTTTTTTTTTTTTTTTGATAAAAAAAA" );
    const auto var1 = std::make_shared< Variant >( referenceSequence, Region( chrom, 106703, 106703 ), "TTT" );
    const auto var2 = std::make_shared< Variant >( referenceSequence, Region( chrom, 106720, 106721 ), "" );

    const variantSet_t variants = {var1, var2};
    const auto padding = 6;

    SetRegions regions;
    regions.insert( Region( chrom, 106701, 106703 ) );
    regions.insert( Region( chrom, 106711, 106713 ) );
    regions.insert( Region( chrom, 106719, 106721 ) );

    const Haplotype theHap( referenceSequence, regions, variants, padding );

    const auto & sequences = theHap.paddedSequences();
    BOOST_REQUIRE_EQUAL( sequences.size(), 2 );
    BOOST_CHECK_EQUAL( "AAAAAACCTTTTTTTTT", sequences[0] );
    BOOST_CHECK_EQUAL( "TTTTTTTTTTTTTTGTAAAAA", sequences[1] );

    const auto & hapRegions = theHap.regions();
    SetRegions expectedRegions;
    expectedRegions.insert( Region( chrom, 106701, 106703 ) );
    expectedRegions.insert( Region( chrom, 106711, 106721 ) );
    BOOST_REQUIRE_EQUAL( hapRegions, expectedRegions );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testVarCalledForSNPs )
{
    using namespace echidna::variant;

    std::string chrom = "1";

    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( chrom, 1000007, 1000016 ), "GCCTCACCC" );
    const varPtr_t var1 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000010, 1000011 ), "A" );
    const varPtr_t var2 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000014, 1000015 ), "G" );

    const variantSet_t variants = {var1};
    const Haplotype haplotype( referenceSequence, referenceSequence->region(), variants, 0 );

    BOOST_CHECK_EQUAL( echidna::caller::Call::VAR, haplotype.getCallType( var1 ) );
    BOOST_CHECK_EQUAL( echidna::caller::Call::REF, haplotype.getCallType( var2 ) );
}

BOOST_AUTO_TEST_CASE( testVarCalledForIndels )
{
    using namespace echidna::variant;

    std::string chrom = "1";

    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( chrom, 1000007, 1000016 ), "GCCTCACCC" );
    const varPtr_t var1 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000010, 1000011 ), "AT" );
    const varPtr_t var2 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000012, 1000014 ), "" );
    const varPtr_t var3 = std::make_shared< Variant >( referenceSequence, Region( chrom, 1000014, 1000015 ), "GC" );

    const variantSet_t variants = {var1, var2};
    const Haplotype haplotype( referenceSequence, referenceSequence->region(), variants, 0 );

    BOOST_CHECK_EQUAL( echidna::caller::Call::VAR, haplotype.getCallType( var1 ) );
    BOOST_CHECK_EQUAL( echidna::caller::Call::VAR, haplotype.getCallType( var2 ) );
    BOOST_CHECK_EQUAL( echidna::caller::Call::REF, haplotype.getCallType( var3 ) );
}

BOOST_AUTO_TEST_CASE( testVarCalledForOverlappingVariants )
{
    using namespace echidna::variant;

    std::string chrom = "1";

    auto refSeq = std::make_shared< ReferenceSequence >( Region( chrom, 1000210, 1000230 ), "ACTTGTTTTAGTTTTGTTTT" );
    const varPtr_t var1 = std::make_shared< Variant >( refSeq, Region( chrom, 1000216, 1000221 ), "" );
    const varPtr_t var2 = std::make_shared< Variant >( refSeq, Region( chrom, 1000219, 1000220 ), "T" );

    const Haplotype haplotype1( refSeq, refSeq->region(), {var1}, 0 );

    BOOST_CHECK_EQUAL( echidna::caller::Call::VAR, haplotype1.getCallType( var1 ) );
    BOOST_CHECK_EQUAL( echidna::caller::Call::UNKNOWN, haplotype1.getCallType( var2 ) );

    const Haplotype haplotype2( refSeq, refSeq->region(), {var2}, 0 );

    BOOST_CHECK_EQUAL( echidna::caller::Call::UNKNOWN, haplotype2.getCallType( var1 ) );
    BOOST_CHECK_EQUAL( echidna::caller::Call::VAR, haplotype2.getCallType( var2 ) );
}