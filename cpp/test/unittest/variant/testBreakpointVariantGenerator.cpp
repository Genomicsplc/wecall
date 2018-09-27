// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <caller/params.hpp>
#include <io/readDataSet.hpp>

#include "variant/breakpointVariantGenerator.hpp"
#include "variant/haplotype.hpp"

using echidna::utils::ReferenceSequence;
using echidna::utils::referenceSequencePtr_t;
using echidna::utils::BasePairSequence;
using echidna::caller::Region;
using echidna::variant::BreakpointVariantGenerator;
using echidna::variant::BreakpointClusterer;
using echidna::variant::BreakpointLocus;
using echidna::variant::breakpointLocusSet_t;
using echidna::variant::Breakpoint;
using echidna::variant::Variant;
using echidna::io::Read;
using echidna::io::ReadDataset;

std::shared_ptr< Read > readFromBreakpoint( const Breakpoint & bp,
                                            const referenceSequencePtr_t & referenceSequence,
                                            const int64_t readLength = 110 )
{
    const int64_t paddingLength = int64_t( readLength ) - int64_t( bp.sequence().size() );
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( echidna::caller::Region( "1", 100, 110 ), "ACGTAAAAGT" );
    if ( bp.isStartBreakpoint() )
    {
        const int64_t start = bp.pos() - paddingLength;
        const auto readSequence =
            echidna::variant::Haplotype::buildHaplotypeSequence(
                referenceSequence, Region( bp.contig(), start, bp.pos() ), bp.getLocalVariants() ) +
            bp.sequence();
        const echidna::utils::QualitySequence qualitySequence( readSequence.size(), 20 );

        const auto cigarString = std::to_string( paddingLength ) + "M" + std::to_string( bp.sequence().size() ) + "S";
        return std::make_shared< Read >( readSequence, qualitySequence, "", cigarString, 0, start, 0, 0, 0, 0, 0,
                                         refSequence );
    }
    else
    {
        const int64_t end = bp.pos() + paddingLength;
        const auto readSequence =
            bp.sequence() + echidna::variant::Haplotype::buildHaplotypeSequence(
                                referenceSequence, Region( bp.contig(), bp.pos(), end ), bp.getLocalVariants() );
        const echidna::utils::QualitySequence qualitySequence( readSequence.size(), 20 );

        const auto cigarString = std::to_string( bp.sequence().size() ) + "S" + std::to_string( paddingLength ) + "M";
        return std::make_shared< Read >( readSequence, qualitySequence, "", cigarString, 0, bp.pos(), 0, 0, 0, 0, 0,
                                         refSequence );
    }
}

BOOST_AUTO_TEST_CASE( testShouldGetZeroClustersForEmptyInput )
{
    const auto paddingDistance = 30;
    const auto clusters = BreakpointClusterer( paddingDistance ).getClusters( {} );
    BOOST_REQUIRE_EQUAL( clusters.size(), 0 );
}

BOOST_AUTO_TEST_CASE( testShouldNotPutBreakpointsIntoClusterIfFarApart )
{
    const auto contig = "20";
    const auto eventStart = 38337632;
    const auto paddingDistance = 30;

    const auto locus1 = std::make_shared< BreakpointLocus >( contig, eventStart, true );
    const auto locus2 = std::make_shared< BreakpointLocus >( contig, eventStart + paddingDistance * 2, true );

    const auto clusters = BreakpointClusterer( paddingDistance ).getClusters( {locus1, locus2} );
    BOOST_REQUIRE_EQUAL( clusters.size(), 2 );
}

BOOST_AUTO_TEST_CASE( testShouldPutBreakpointsIntoClusterIfWithin2TimesDistance )
{
    const auto contig = "20";
    const auto eventStart = 38337632;
    const auto paddingDistance = 30;

    const auto locus1 = std::make_shared< BreakpointLocus >( contig, eventStart, true );
    const auto locus2 = std::make_shared< BreakpointLocus >( contig, eventStart + paddingDistance * 2 - 1, true );

    const auto clusters = BreakpointClusterer( paddingDistance ).getClusters( {locus1, locus2} );
    BOOST_REQUIRE_EQUAL( clusters.size(), 1 );
}

BOOST_AUTO_TEST_CASE( testShouldClusterTogetherIfPositionContainedInMateRegionOfOther )
{
    const auto contig = "20";
    const auto eventStart = 38337632;
    const auto eventEnd = eventStart + 100000;
    const auto paddingDistance = 30;

    const auto locus1 = std::make_shared< BreakpointLocus >( contig, eventStart, true );
    locus1->addMateRegion( Region( contig, eventEnd + paddingDistance - 1, eventEnd + paddingDistance * 2 ), 0 );

    const auto locus2 = std::make_shared< BreakpointLocus >( contig, eventEnd, true );

    const auto clusters = BreakpointClusterer( paddingDistance ).getClusters( {locus1, locus2} );
    BOOST_REQUIRE_EQUAL( clusters.size(), 1 );
}

BOOST_AUTO_TEST_CASE( testShouldNotClusterTogetherIfPositionBeforePaddedMateRegionOfOther )
{
    const auto contig = "20";
    const auto eventStart = 38337632;
    const auto eventEnd = eventStart + 100000;
    const auto paddingDistance = 30;

    const auto locus1 = std::make_shared< BreakpointLocus >( contig, eventStart, true );
    locus1->addMateRegion( Region( contig, eventEnd + paddingDistance, eventEnd + paddingDistance * 2 ), 0 );

    const auto locus2 = std::make_shared< BreakpointLocus >( contig, eventEnd, true );

    const auto clusters = BreakpointClusterer( paddingDistance ).getClusters( {locus1, locus2} );
    BOOST_REQUIRE_EQUAL( clusters.size(), 2 );
}

BOOST_AUTO_TEST_CASE( testBreakpointClusteringShouldClusterTogether )
{
    const auto contig = "20";
    const auto eventStart = 38337632;
    const auto eventEnd = 38337769;

    const auto locus1 = std::make_shared< BreakpointLocus >( contig, eventStart + 15, true );
    locus1->addMateRegion( Region( contig, 38336681, 38338062 ) );

    const auto locus2 = std::make_shared< BreakpointLocus >( contig, eventEnd - 40, false );
    locus2->addMateRegion( Region( contig, 38337506, 38337988 ) );

    const auto locus3 = std::make_shared< BreakpointLocus >( contig, eventEnd, false );
    locus3->addMateRegion( Region( contig, 38337287, 38338046 ) );

    BOOST_CHECK_EQUAL( 1, BreakpointClusterer( 30 ).getClusters( {locus1, locus2, locus3} ).size() );
}

BOOST_AUTO_TEST_CASE( shouldClusterStartAndEndBreakpointsTogetherEvenIfThereIsABreakpointInbetween )
{
    const auto contig = "3";

    const auto locus1 = std::make_shared< BreakpointLocus >( contig, 190290650, true );
    locus1->addMateRegion( Region( contig, 190290036, 190291134 ) );

    const auto locus2 = std::make_shared< BreakpointLocus >( contig, 190290670, true );
    locus2->addMateRegion( Region( contig, 190290198, 190290976 ) );

    const auto locus3 = std::make_shared< BreakpointLocus >( contig, 190290858, true );

    const auto locus4 = std::make_shared< BreakpointLocus >( contig, 190290936, false );
    locus4->addMateRegion( Region( contig, 190290407, 190290958 ) );

    const auto clusters = BreakpointClusterer( 30 ).getClusters( {locus1, locus2, locus3, locus4} );
    BOOST_CHECK_EQUAL( 1, clusters.size() );
}

BOOST_AUTO_TEST_CASE( testShouldFindDeletionAt373682 )
{
    const auto contig = "20";
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( contig, 373540, 373900 ),
        "GAAGGGCAGAGAGACTCGCCTCCTGTTCCCTCTCCAGATTCCTGGGGGCAGTCAAGATGTGTCAGGGAGTGCACTAAGCTGCCAGTTACACAGGAGTTTTCTGTGGAGAA"
        "AGGAGTGTGACCCCATGGCATTTTAAAAAACTTTTTATCTTGAAATAATTTTAGACTTTTAGAAAACCTACAAAAATAGTTCAAAGAGTTTCTGCATATCCTTTAACCAG"
        "TGCTCTCCAATGTTAACACCTGACGTAGCTATGGTACAATTACCCAAACTATTAACTAAGCCACAGATTGATTCCCACTTCCCGAGTTTCCCACTAACACCCCTTGCTGT"
        "GCAGGATCCAGTGAGGATCCCAATTTAGTC" );

    BreakpointVariantGenerator breakpointVariantGenerator( referenceSequence,
                                                           echidna::caller::params::defaults::defaultBreakpointKmerSize,
                                                           echidna::caller::params::defaults::maxBreakpointKmerSize );

    const BasePairSequence added = "AAGGGGTGTTAGTGGGTACTAAGAAGTCAACCTTGGTAGAGT";
    const auto eventStart = 373682;
    const auto eventEnd = 373807;

    const auto startBreakpoint = std::make_shared< Breakpoint >(
        contig, eventStart, true,
        added + referenceSequence->subseq( Region( contig, eventEnd, eventEnd + 24 ) ).sequence() );

    const auto endBreakpoint = std::make_shared< Breakpoint >(
        contig, eventEnd, false,
        referenceSequence->subseq( Region( contig, eventStart - 32, eventStart ) ).sequence() + added );

    auto startBpLocus = std::make_shared< BreakpointLocus >( contig, eventStart, true );
    auto endBpLocus = std::make_shared< BreakpointLocus >( contig, eventEnd, false );

    ReadDataset readDataset( {""}, referenceSequence->region() );
    readDataset.insertRead( "", readFromBreakpoint( *startBreakpoint, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *endBreakpoint, referenceSequence ) );
    const auto allReads = readDataset.getAllReads( 0 );

    const auto variants = breakpointVariantGenerator.getVariantCandidates( {startBpLocus, endBpLocus}, allReads );
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    const std::vector< echidna::variant::varPtr_t > vecVariant( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( *vecVariant[0], Variant( referenceSequence, Region( contig, eventStart, eventEnd ), added ) );

    for ( const auto & variant : variants )
    {
        BOOST_CHECK( variant->isFromBreakpoint() );
    }
}

BOOST_AUTO_TEST_CASE( testWithLongOverlappingSoftClippedSequenceAt135116 )
{
    const auto contig = "20";
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( contig, 135000, 135300 ),
        "TCATAATGAATAATTCTGTATTTTGATTTGTTGATGGTTAAATGAATCTAAGAAGGATCTGATAATCTAGTATTAGACGCATAGTACAGTAGTGAACTTTCAGCCTCTAG"
        "AGAAGTGCTCTCCAATAGAAGTAACTAAAAAGATGGAAATGTTCTAAATCTGTACTATACAATATGTTAGCCACTGGTCATGTGGTTGTAAAGCACTTGATATGTGCTCT"
        "ATAATATTTTGAGCAAGTTTCTTATGTTTTCTGTGCCTTTGAGTTTCTCATTTGTAGAATGGAGATAATAATATCTACCT" );

    BreakpointVariantGenerator breakpointVariantGenerator( referenceSequence,
                                                           echidna::caller::params::defaults::defaultBreakpointKmerSize,
                                                           echidna::caller::params::defaults::maxBreakpointKmerSize );
    const BasePairSequence added1 = "TTGCATTCACAGGAAGAATATAAATAACT";
    const BasePairSequence added2 = "TTGCATTCACAGGAAGACTATAAATAACT";
    const auto eventStart = 135116;
    const auto eventEnd = 135120;

    const auto startBreakpoint = std::make_shared< Breakpoint >(
        contig, eventStart, true, added1 + BasePairSequence( "TCCAATAGAAGTAACTAAAAAGATGGAAAT" ) );
    const auto startBreakpoint2 = std::make_shared< Breakpoint >(
        contig, eventStart, true, added2 + BasePairSequence( "TCCAATAGAAGTAACTAAAAAGATGGAAATGTTCT" ) );
    const auto endBreakpoint = std::make_shared< Breakpoint >(
        contig, eventEnd, false, BasePairSequence( "ATAGTACAGTAGTGAACTTTCAGCCTCTAGAGAAGT" ) + added2 );

    auto startBpLocus = std::make_shared< BreakpointLocus >( contig, eventStart, true );
    auto endBpLocus = std::make_shared< BreakpointLocus >( contig, eventEnd, false );

    ReadDataset readDataset( {""}, referenceSequence->region() );
    readDataset.insertRead( "", readFromBreakpoint( *startBreakpoint, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *startBreakpoint2, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *endBreakpoint, referenceSequence ) );
    const auto allReads = readDataset.getAllReads( 0 );

    const auto variants = breakpointVariantGenerator.getVariantCandidates( {startBpLocus, endBpLocus}, allReads );
    BOOST_REQUIRE_EQUAL( variants.size(), 2 );
    const std::vector< echidna::variant::varPtr_t > vecVariant( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( *vecVariant[0], Variant( referenceSequence, Region( contig, eventStart, eventEnd ), added1 ) );

    BOOST_CHECK_EQUAL( *vecVariant[1], Variant( referenceSequence, Region( contig, eventStart, eventEnd ), added2 ) );

    for ( const auto & variant : variants )
    {
        BOOST_CHECK( variant->isFromBreakpoint() );
    }
}

BOOST_AUTO_TEST_CASE( testMediumInsertionAt1745352 )
{
    const auto contig = "20";
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( contig, 1745200, 1745450 ),
        "TCTCCCTGTGCCTCAGGCTAAGTAGGGAGGCCTGAGGAGAGGGAGCAAGACGGGGGAACGGCCAGTAGGTGGAGCAGTCAGAACACACAAGACATTTTCCATCAAGTGTG"
        "CTGTTTCATAAGGACACGGTTCGCGGAGCCACAAGACCATTACACTAATCACAGATCCTCAGATCAGATATAATAATTATTAAAATTTGAAATATTGCAAGAATTACCAA"
        "AATGTGACAGAGACACGGAGTGGGCACATG" );

    BreakpointVariantGenerator breakpointVariantGenerator( referenceSequence,
                                                           echidna::caller::params::defaults::defaultBreakpointKmerSize,
                                                           echidna::caller::params::defaults::maxBreakpointKmerSize );
    const BasePairSequence tailAdded = "TCAAAGAT";
    const auto eventStart = 1745352;
    const auto eventEnd = 1745352;
    const BasePairSequence added1 = BasePairSequence( "CAATAATCACA" ) + tailAdded;
    const BasePairSequence added2 =
        referenceSequence->subseq( Region( contig, eventStart, eventStart + 11 ) ).sequence() + tailAdded;

    const auto endBreakpoint = std::make_shared< Breakpoint >(
        contig, eventEnd, false,
        referenceSequence->subseq( Region( contig, 1745332, eventStart ) ).sequence() + added1 );
    const auto startBreakpoint = std::make_shared< Breakpoint >(
        contig, eventStart + 11, true,
        tailAdded + referenceSequence->subseq( Region( contig, eventEnd, 1745377 ) ).sequence() );

    auto startBpLocus = std::make_shared< BreakpointLocus >( contig, eventStart + 11, true );
    auto endBpLocus = std::make_shared< BreakpointLocus >( contig, eventEnd, false );

    ReadDataset readDataset( {""}, referenceSequence->region() );
    readDataset.insertRead( "", readFromBreakpoint( *startBreakpoint, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *endBreakpoint, referenceSequence ) );
    const auto allReads = readDataset.getAllReads( 0 );

    const auto variants = breakpointVariantGenerator.getVariantCandidates( {startBpLocus, endBpLocus}, allReads );
    BOOST_REQUIRE_EQUAL( variants.size(), 2 );
    const std::vector< echidna::variant::varPtr_t > vecVariant( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( *vecVariant[0], Variant( referenceSequence, Region( contig, eventStart, eventEnd ), added1 ) );
    BOOST_CHECK_EQUAL( *vecVariant[1], Variant( referenceSequence, Region( contig, eventStart, eventEnd ), added2 ) );

    for ( const auto & variant : variants )
    {
        BOOST_CHECK( variant->isFromBreakpoint() );
    }
}

BOOST_AUTO_TEST_CASE( testBobbleWithEndBreakpointBeforeStartBreakpointAt1778684 )
{
    auto contig = "20";
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( contig, 1778600, 1778800 ),
        "CTACTTAGGAGGATGAGGCAGGAGAATCACTTGCACCCGGGAGGTGGAGGCTGCAGTGAGCTGAGATCCCACCATTGCCCTCCATTGGAGGCTGAGGTGGAAGGATTGCT"
        "TGAATTCAGGTGTTCAAGACCAGCCTGGACAACATAGCAAGACCCCGCCTCTAAAAAAAAATTAGCTGGGCATTGTAGAATACCTGTAGT" );

    BreakpointVariantGenerator breakpointVariantGenerator( referenceSequence,
                                                           echidna::caller::params::defaults::defaultBreakpointKmerSize,
                                                           echidna::caller::params::defaults::maxBreakpointKmerSize );
    const auto eventStart = 1778684;
    const auto eventEnd = 1778684;
    const BasePairSequence added =
        "TTGCCTAGGTGACAGAGTGAGGCTCCATCTTAAGAAAAAAAAAAGCCTGTTAAGAGGATGAAAAGAGCAGCCAGGCATTGACTGTGGTGCCTCACGCCTATAATCCCAGT"
        "GCTT";

    const auto endBreakpoint =
        std::make_shared< Breakpoint >( contig, eventEnd, false, added.substr( added.size() - 75, 75 ) );
    const auto startBreakpoint = std::make_shared< Breakpoint >( contig, eventStart + 3, true, added.substr( 3, 80 ) );

    auto startBpLocus = std::make_shared< BreakpointLocus >( contig, eventStart + 3, true );
    auto endBpLocus = std::make_shared< BreakpointLocus >( contig, eventEnd, false );

    ReadDataset readDataset( {""}, referenceSequence->region() );
    readDataset.insertRead( "", readFromBreakpoint( *startBreakpoint, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *endBreakpoint, referenceSequence ) );
    const auto allReads = readDataset.getAllReads( 0 );

    const auto variants = breakpointVariantGenerator.getVariantCandidates( {startBpLocus, endBpLocus}, allReads );
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    const std::vector< echidna::variant::varPtr_t > vecVariant( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( *vecVariant[0], Variant( referenceSequence, Region( contig, eventStart, eventEnd ), added ) );

    for ( const auto & variant : variants )
    {
        BOOST_CHECK( variant->isFromBreakpoint() );
    }
}

BOOST_AUTO_TEST_CASE( testBobbleWithEndBreakpointBeforeStartBreakpointAt5291304 )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( "20", 5291100, 5291500 ),
        "ACTGTGCCCAGTGATAAAAGCTCCCTTCTCCAAAATCAGCATAATTGGGGCAAAGGTGAGAGTTTGGCCATCCATCTTTCTTAGCATTTTGCTTTATTTTTTTCAATTAG"
        "TTTAGAACCTGATCAGAGAAGAGGATCCTTACTCTTTCTCCCTATCCTACCACCTGCCACCTGTACAGCATTTCCAGAGCAGCCCTCATCAATGACCTCCCACCTCTGGA"
        "AGAACTCTCTGTAGAGCTGTGATCTATGCTGGCTCCACTAGTCCCTCAGCAGGAATATGCTCTAGTTGCCACAGTGAAAGCTTGCTTAATTACTTCTTGGGTGTTACCTG"
        "CACATAAATCTTAGTCTCAGAGTCTGATTTGGGAGTAACCCAAGCTAAGACTCAGATATAAGTCCCTGTG" );

    BreakpointVariantGenerator breakpointVariantGenerator( referenceSequence,
                                                           echidna::caller::params::defaults::defaultBreakpointKmerSize,
                                                           echidna::caller::params::defaults::maxBreakpointKmerSize );

    const BasePairSequence altSeq = "TTATGGGGAGCTGGGGTATAAATACCCCAGCTCCCTCCT";
    const auto eventStart = 5291304;
    const auto eventEnd = 5291304;
    const auto breakpoint1 = std::make_shared< Breakpoint >(
        "20", eventEnd, false, referenceSequence->subseq( Region( "20", 5291200, eventStart - 5 ) ).sequence() +
                                   BasePairSequence( "TAATGACT" ) + altSeq );

    const auto breakpoint2 = std::make_shared< Breakpoint >(
        "20", eventStart + 2, true, altSeq + referenceSequence->subseq( Region( "20", 5291304, 5291400 ) ).sequence() );

    auto startBpLocus = std::make_shared< BreakpointLocus >( "20", eventStart + 2, true );
    auto endBpLocus = std::make_shared< BreakpointLocus >( "20", eventEnd, false );

    ReadDataset readDataset( {""}, referenceSequence->region() );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint1, referenceSequence, 150 ) );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint2, referenceSequence, 150 ) );
    const auto allReads = readDataset.getAllReads( 0 );

    const auto variants = breakpointVariantGenerator.getVariantCandidates( {startBpLocus, endBpLocus}, allReads );

    BOOST_REQUIRE( variants.size() == 2 );
    const std::vector< echidna::variant::varPtr_t > vecVariant( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( *vecVariant[0], Variant( referenceSequence, Region( "20", eventStart - 5, eventEnd ),
                                                BasePairSequence( "TAATGACT" ) + altSeq ) );

    BOOST_CHECK_EQUAL( *vecVariant[1], Variant( referenceSequence, Region( "20", eventStart, eventEnd ),
                                                BasePairSequence( "AC" ) + altSeq ) );

    for ( const auto & variant : variants )
    {
        BOOST_CHECK( variant->isFromBreakpoint() );
    }
}

BOOST_AUTO_TEST_CASE( testBobbleWithEndBreakpointBeforeStartBreakpointAt6645839 )
{
    const auto contig = "20";
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( contig, 6645700, 6645900 ),
        "ATTTCCAAGAGTTAGTATTGAAAATGTAAATTATCTTATTATTTCATATTGATTACATGTTGAAGTAGTATTATTTTGGATATATTGAGTTAAATAAAATATATTGTTAA"
        "AATTCATTTTCTCTATTTTTTTTTTTTTTTTGCTTATTTGTATTTCTGCTGTTGCTCTAGATTAACCAAACTGATCTCCACATTCTCCAT" );
    const auto eventStart = 6645839;
    const auto eventEnd = 6645839;

    BreakpointVariantGenerator breakpointVariantGenerator( referenceSequence,
                                                           echidna::caller::params::defaults::defaultBreakpointKmerSize,
                                                           echidna::caller::params::defaults::maxBreakpointKmerSize );

    const BasePairSequence added = "TTTTACTTTTTTAATGTGGCTACTAGAAAACTTGAAATTCCGTATGTGGCTTGCATTGCATTTCTA";
    const auto breakpoint1 = std::make_shared< Breakpoint >(
        contig, eventEnd, false,
        referenceSequence->subseq( Region( contig, eventStart - 12, eventStart ) ).sequence() + added );

    // This breakpoint has some junk on the end.
    const auto breakpoint2 = std::make_shared< Breakpoint >(
        contig, eventStart + 2, true, added.substr( 2, 61 ) + BasePairSequence( "TTTTTTCTTTTTTGT" ) );

    auto startBpLocus = std::make_shared< BreakpointLocus >( contig, eventStart + 2, true );
    auto endBpLocus = std::make_shared< BreakpointLocus >( contig, eventStart - 12, false );

    ReadDataset readDataset( {""}, referenceSequence->region() );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint1, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint2, referenceSequence ) );
    const auto allReads = readDataset.getAllReads( 0 );

    const auto variants = breakpointVariantGenerator.getVariantCandidates( {startBpLocus, endBpLocus}, allReads );

    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    const std::vector< echidna::variant::varPtr_t > vecVariant( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( *vecVariant[0], Variant( referenceSequence, Region( contig, eventStart, eventEnd ), added ) );

    for ( const auto & variant : variants )
    {
        BOOST_CHECK( variant->isFromBreakpoint() );
    }
}

BOOST_AUTO_TEST_CASE( testFindsDeletionAt9696769 )
{
    // 1000G known structural variant: DEL_pindel_49611
    const auto contig = "20";
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( contig, 9696620, 9696960 ),
        "CAACATGGCACATGTATGCATATGTAACAAACCTGCACATTGTGCAGATGTACCCTAGAACTTAAAGTGTAATAAAAAATATGTACATAAAAATCAAAATCAAAGAAAGA"
        "ACATGCAGTAGCTGAAAAAAAATATCTTCTCAAAAGCATACCCTAAAACTTAGAGTATAATAAAAAAAAAAAAATTAAAAAAAAAAAAAGCATGCTCTATGTTTTAAACT"
        "ATTATTGCTAGGATCACTAGGACTTAGTAAAAAGCAATGCCTTACACAGGCAACAATAAATTCTAGAAGGACAGTTCTGACAGCTGCTTTTCAATTAGTTCAAAATAGGC"
        "TAAAATATAA" );

    BreakpointVariantGenerator breakpointVariantGenerator( referenceSequence,
                                                           echidna::caller::params::defaults::defaultBreakpointKmerSize,
                                                           echidna::caller::params::defaults::maxBreakpointKmerSize );
    const auto eventStart = 9696761;
    const auto eventEnd = 9696815;
    const auto overlap = 8;

    const auto breakpoint1 = std::make_shared< Breakpoint >(
        contig, eventStart + overlap, true,
        referenceSequence->subseq( Region( contig, eventEnd + overlap, 9696850 ) ).sequence() +
            BasePairSequence( "GGATCACTAGGACTTAGTA" ) );

    const auto breakpoint2 = std::make_shared< Breakpoint >(
        contig, eventEnd, false, BasePairSequence( "AAAATCAAAGAAAGA" ) +
                                     referenceSequence->subseq( Region( contig, 9696730, eventStart ) ).sequence() );

    auto startBpLocus = std::make_shared< BreakpointLocus >( contig, eventStart + overlap, true );
    auto endBpLocus = std::make_shared< BreakpointLocus >( contig, eventEnd, false );

    ReadDataset readDataset( {""}, referenceSequence->region() );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint1, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint2, referenceSequence ) );
    const auto allReads = readDataset.getAllReads( 0 );

    const auto variants = breakpointVariantGenerator.getVariantCandidates( {startBpLocus, endBpLocus}, allReads );
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    const std::vector< echidna::variant::varPtr_t > vecVariant( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( *vecVariant[0], Variant( referenceSequence, Region( contig, eventStart, eventEnd ), "" ) );

    for ( const auto & variant : variants )
    {
        BOOST_CHECK( variant->isFromBreakpoint() );
    }
}

BOOST_AUTO_TEST_CASE( testFindsComplexDeletionAt11398595 )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( "20", 11398400, 11398700 ),
        "CAGACCTATAATATTTTTCTTCTAACTGAAGAATTTTTTAGCATTTCTTGCAAGTCAGGCCTGCTAGTGATTAATTCTCTCAGTATTTTTTGTCTGATAAAGTCAATATC"
        "CCTCATTTTTGAAACATAGCTTTACTAGATATAGAGTTCCAGATTAATGTTTGCTTTTTTTAATCACTTT"
        "AATTAATTTATTTTAAATTAATTAAATTTAAAATTTAATTAATTTTAATTCTTGCTTTCA"
        "TGGTTTCTGGTGAGAAGTCTGCTGTAATTTTAACCTTTTTCTCTATAGATAAGGTGATTT" );

    BreakpointVariantGenerator breakpointVariantGenerator( referenceSequence,
                                                           echidna::caller::params::defaults::defaultBreakpointKmerSize,
                                                           echidna::caller::params::defaults::maxBreakpointKmerSize );

    const auto breakpoint1 =
        std::make_shared< Breakpoint >( "20", 11398595, true, "CTGTCTTCTTGCTTTCATGGTTTCTGGTGAGAAGTCTGCTGTAATTTTAA" );

    const auto breakpoint2 =
        std::make_shared< Breakpoint >( "20", 11398628, false, "ATTAATGTTTGCTTTTTTTAATCACTTTAATTAATTTATTTTACTGTC" );

    auto startBpLocus = std::make_shared< BreakpointLocus >( "20", 11398595, true );
    auto endBpLocus = std::make_shared< BreakpointLocus >( "20", 11398628, false );

    ReadDataset readDataset( {""}, referenceSequence->region() );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint1, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint2, referenceSequence ) );
    const auto allReads = readDataset.getAllReads( 0 );

    const auto variants = breakpointVariantGenerator.getVariantCandidates( {startBpLocus, endBpLocus}, allReads );
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    const std::vector< echidna::variant::varPtr_t > vecVariant( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( *vecVariant[0], Variant( referenceSequence, Region( "20", 11398595, 11398628 ), "CTGTC" ) );
}

BOOST_AUTO_TEST_CASE( testShouldNotFindVariantsAsBelievedToBeGenuineBreakpoints )
{
    const auto referenceSequence =
        std::make_shared< ReferenceSequence >( Region( "20", 12359661, 12359771 ),
                                               "TACGGAAGAATCTATAATTATATTAGGATTTTACAGGGAGAATATTATCTCCCTTCCTAAAGCACCTTCCT"
                                               "AAAGGAGGAGGTATCTCAATGTGCAAAGTGAACCTAATT" );

    BreakpointVariantGenerator breakpointVariantGenerator( referenceSequence,
                                                           echidna::caller::params::defaults::defaultBreakpointKmerSize,
                                                           echidna::caller::params::defaults::maxBreakpointKmerSize );

    const auto breakpoint1 = std::make_shared< Breakpoint >(
        "20", 12359712, true, "TTAGCTCTATAACTATATTTAGAGACTAAAATACCTGGTAATATCCCGTTCAGTCGCCCCTATTGATTTAGGAATTCAT" );

    const auto breakpoint2 = std::make_shared< Breakpoint >(
        "20", 12359721, false, "GCAACTATTTATTCTACCACTCACCAGAAAAGATTGCTTACCATAGCAACCCCTGCATTAGCATTCAGCTTCTTGGTTACAG" );

    auto startBpLocus = std::make_shared< BreakpointLocus >( "20", 12359712, true );
    auto endBpLocus = std::make_shared< BreakpointLocus >( "20", 12359721, false );

    ReadDataset readDataset( {""}, referenceSequence->region() );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint1, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint2, referenceSequence ) );
    const auto allReads = readDataset.getAllReads( 0 );

    const auto variants = breakpointVariantGenerator.getVariantCandidates( {startBpLocus, endBpLocus}, allReads );
    BOOST_REQUIRE_EQUAL( variants.size(), 0 );

    for ( const auto & variant : variants )
    {
        BOOST_CHECK( variant->isFromBreakpoint() );
    }
}

BOOST_AUTO_TEST_CASE( testShouldFind1000GTrueDeletionAt6150945 )
{
    // 1000G known structural variant: DEL_pindel_49512
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( "20", 6150790, 6151350 ),
        "ATCTTTCTGAGCTGCTGGGATAGTACATCTCATTCATTGTTGGAGTTTCCTGTCTGCCCCCATCCCAAGTCAGGTATAAGCCTCAGAGCAGACATTTGCAGAGAATCCCA"
        "TGGTATCCATTAACAAGCTCATGGTTTAAGTGCACATTATCTTTTAAGAGTGGACATCCGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAG"
        "GCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCCGGCTAAAACGGTGAAACCCCGTCTCTACTGAAAATACAAAAAATTAGCCGGGCGTAGTGGCGGGCGCCTGTA"
        "GTCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCCCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAG"
        "ACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAGTGGACATCCTTTCTGTGGAAATGAGAATCTTTGTGGAACCCTGGTGCTTTGGTCCAGT"
        "TGTACTGGGG" );

    BreakpointVariantGenerator breakpointVariantGenerator( referenceSequence,
                                                           echidna::caller::params::defaults::defaultBreakpointKmerSize,
                                                           echidna::caller::params::defaults::maxBreakpointKmerSize );

    const auto breakpoint1 =
        std::make_shared< Breakpoint >( "20", 6150959, true, "TTTCTGTGGAAATGAGAATCTTTGTGGAACCCTGGTGCT" );

    const auto breakpoint2 =
        std::make_shared< Breakpoint >( "20", 6151277, false, "ATCCATTAACAAGCTCATGGTTTAAGTGCACATTATCTTTT" );

    auto startBpLocus = std::make_shared< BreakpointLocus >( "20", 6150959, true );
    auto endBpLocus = std::make_shared< BreakpointLocus >( "20", 6151277, false );

    ReadDataset readDataset( {""}, referenceSequence->region() );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint1, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint2, referenceSequence ) );
    const auto allReads = readDataset.getAllReads( 0 );

    const auto variants = breakpointVariantGenerator.getVariantCandidates( {startBpLocus, endBpLocus}, allReads );
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    const std::vector< echidna::variant::varPtr_t > vecVariant( variants.begin(), variants.end() );
    BOOST_CHECK_EQUAL( *vecVariant[0], Variant( referenceSequence, Region( "20", 6150945, 6151277 ), "" ) );

    for ( const auto & variant : variants )
    {
        BOOST_CHECK( variant->isFromBreakpoint() );
    }
}

BOOST_AUTO_TEST_CASE( testShouldFind1000GTrueDeletionAt19170061 )
{
    // 1000G known structural variant: DEL_pindel_49995
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( "20", 19170000, 19170450 ),
        "AACAAGAAATATATCTCCATTTATTTATACTTCTTTGACCTTTTCAGAGTTTTGTAGTTTTCCTCATAAATTTCTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCG"
        "CTCTGTCGCCCAGGCTGGAGTGCAGTGGCGCGATCTCGGCTCACTGCAAGCTCCGCCTCCTGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACA"
        "GGCGCCCGCCACTACGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCATGATCCGCCCGCCTCG"
        "GCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCCCCTCATAAATTTCTTATACATATTTTGTTAAGTTTTTATGTAAGTATTTTATTTTTTAACACT"
        "AATGTAAATA" );

    BreakpointVariantGenerator breakpointVariantGenerator( referenceSequence,
                                                           echidna::caller::params::defaults::defaultBreakpointKmerSize,
                                                           echidna::caller::params::defaults::maxBreakpointKmerSize );

    const auto breakpoint1 =
        std::make_shared< Breakpoint >( "20", 19170076, true, "ATACATATTTTGTTAAGTTTTTATGTAAGTATTTTATTTTT" );

    const auto breakpoint2 =
        std::make_shared< Breakpoint >( "20", 19170377, false, "TATTTATACTTCTTTGACCTTTTCAGAGTTTTGTAGTTTT" );

    auto startBpLocus = std::make_shared< BreakpointLocus >( "20", 19170076, true );
    auto endBpLocus = std::make_shared< BreakpointLocus >( "20", 19170377, false );

    ReadDataset readDataset( {""}, referenceSequence->region() );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint1, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint2, referenceSequence ) );
    const auto allReads = readDataset.getAllReads( 0 );

    const auto variants = breakpointVariantGenerator.getVariantCandidates( {startBpLocus, endBpLocus}, allReads );
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    const std::vector< echidna::variant::varPtr_t > vecVariant( variants.begin(), variants.end() );
    BOOST_CHECK_EQUAL( *vecVariant[0], Variant( referenceSequence, Region( "20", 19170061, 19170377 ), "" ) );

    for ( const auto & variant : variants )
    {
        BOOST_CHECK( variant->isFromBreakpoint() );
    }
}

BOOST_AUTO_TEST_CASE( testShouldFindDeletionAt38337632 )
{
    // 1000G known structural variant: DEL_pindel_50373
    const auto contig = "20";
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( contig, 38337490, 38337930 ),
        "AAGCCATTTCTTCCTTGTAAGTTTCCACTGAGGGATTGTACTCGCCAAAGAAGTACAAATGCATTTTAGATCTACAATCATTTATCAATTATCAAACACCTATCAAAATC"
        "ATAGGGTTGGCTCTATATTATGTTGGAAGAATAAAAATATGACTCGAGAACCCGGGAAGGCGGAGCTTGCAGTGAGCCGAGATCCCGCCACTGCACTCCAGCCTGGGCGA"
        "CAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATATGACTCGAACAGGCCCTTGACCTTGGTGGACTTTAAGTGAGGGA"
        "GAAATAAATTATCCCAGATATTTGCAATGCTAGAAAGAGATGCAAAGCTAAGATGGGTAAAGAAGAAGTAGTAATTAACACCTTGAAGACTATTTAGCAGTTCTCCAAG"
        "A" );

    BreakpointVariantGenerator breakpointVariantGenerator( referenceSequence,
                                                           echidna::caller::params::defaults::defaultBreakpointKmerSize,
                                                           echidna::caller::params::defaults::maxBreakpointKmerSize );
    const auto eventStart = 38337632;
    const auto eventEnd = 38337769;

    const auto breakpoint1 =
        std::make_shared< Breakpoint >( contig, eventStart + 15, true, "ACAGGCCCTTGACCTTGGTGGACTTTAAGTGAGGGAGAAATA" );

    const auto breakpoint2 = std::make_shared< Breakpoint >( contig, eventEnd - 40, false, "ATGCCATGAAAAAAAAAAA" );

    const auto breakpoint3 =
        std::make_shared< Breakpoint >( contig, eventEnd, false, "TATCAAAATCATAGGGTTGGCTCTATATTATGTTGGAAGAAT" );

    auto startBpLocus = std::make_shared< BreakpointLocus >( contig, eventStart + 15, true );
    auto endBpLocus = std::make_shared< BreakpointLocus >( contig, eventEnd - 40, false );
    auto endBpLocus2 = std::make_shared< BreakpointLocus >( contig, eventEnd, false );

    ReadDataset readDataset( {""}, referenceSequence->region() );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint1, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint2, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint3, referenceSequence ) );
    const auto allReads = readDataset.getAllReads( 0 );

    const auto variants =
        breakpointVariantGenerator.getVariantCandidates( {startBpLocus, endBpLocus, endBpLocus2}, allReads );
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    const std::vector< echidna::variant::varPtr_t > vecVariant( variants.begin(), variants.end() );
    BOOST_CHECK_EQUAL( *vecVariant[0], Variant( referenceSequence, Region( contig, eventStart, eventEnd ), "" ) );
}

BOOST_AUTO_TEST_CASE( testShouldFindDeletionAt48543735 )
{
    // 1000G known structural variant: DEL_pindel_50691
    const auto contig = "20";
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( contig, 48543590, 48543900 ),
        "TAAGGAGTTCGAGACCAGCCTGGCCAACATGGTGAAACCCTGTCTCTACTAAAAATACAAAAAAATTAGCTGGACGTGGTGGCACACACTTGAAATCCCAGCTACTAGGG"
        "AGGCTGAGGGAGAAGAATTACTTGAACCGGGAGGCAAAGGTTGCAGTGAGCCGAGATGGTACCACTGCACTCCAGCCTGGGAACAAAGTGAGACTCCATCTCAAAAAAAA"
        "AAATAAATAAATAAATAACAGAAGATCTGGACAGACACTTCACAGAGGTAAAACATACAAAATGGCCAATAAACAAATGAAAAGGTACTC" );

    BreakpointVariantGenerator breakpointVariantGenerator( referenceSequence,
                                                           echidna::caller::params::defaults::defaultBreakpointKmerSize,
                                                           echidna::caller::params::defaults::maxBreakpointKmerSize );
    const auto eventStart = 48543735;
    const auto eventEnd = 48543815;

    const auto breakpoint1 = std::make_shared< Breakpoint >(
        contig, eventStart + 2, true,
        referenceSequence->subseq( Region( contig, eventEnd + 2, eventEnd + 16 ) ).sequence() +
            BasePairSequence( "G" ) +
            referenceSequence->subseq( Region( contig, eventEnd + 17, 48543864 ) ).sequence() );

    const auto breakpoint2 = std::make_shared< Breakpoint >(
        contig, eventEnd, false, referenceSequence->subseq( Region( contig, 48543689, eventStart ) ).sequence() );

    const auto nearSnp1 =
        std::make_shared< Variant >( referenceSequence, Region( contig, eventEnd + 16, eventEnd + 17 ), "G" );
    breakpoint2->addLocalVariants( {nearSnp1} );

    const auto breakpoint3 = std::make_shared< Breakpoint >(
        contig, eventEnd, false, referenceSequence->subseq( Region( contig, 48543689, eventStart ) ).sequence() );

    auto startBpLocus = std::make_shared< BreakpointLocus >( contig, eventStart + 2, true );
    auto endBpLocus = std::make_shared< BreakpointLocus >( contig, eventEnd, false );
    startBpLocus->add( breakpoint1 );
    endBpLocus->add( breakpoint2 );
    endBpLocus->add( breakpoint3 );

    ReadDataset readDataset( {""}, referenceSequence->region() );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint1, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint2, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint3, referenceSequence ) );
    const auto allReads = readDataset.getAllReads( 0 );

    const auto variants = breakpointVariantGenerator.getVariantCandidates( {startBpLocus, endBpLocus}, allReads );
    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    const std::vector< echidna::variant::varPtr_t > vecVariant( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( *vecVariant[0], Variant( referenceSequence, Region( contig, eventStart, eventEnd ), "" ) );
}

BOOST_AUTO_TEST_CASE( testShouldFindDeletionAtChr22_50664005 )
{
    // 1000G known structural variant: DEL_pindel_52851
    const auto contig = "22";
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( contig, 50663900, 50664150 ),
        "GAGTTTCCCACGCCAGCAAAGCCTGGGCTTCAGGCTGAGCCACCCAGACCCAGGAGGAGAAGTCAGAGCCCCTGCCTGGCATCTGGCTCCGGACATTCATCTCAGCTCCC"
        "CCTTGGCACCCACCCTTCATCTCAACTCCCCGTTGACACCCACCCTTCATCTCAACTCCCCCTTGGCACCCACCCAGGCCACCGCCAGGACAGAAAATGCCAGAAGGGAG"
        "GCCCCCATCTTGCACAGCCCTCATGTCCGC" );

    BreakpointVariantGenerator breakpointVariantGenerator( referenceSequence,
                                                           echidna::caller::params::defaults::defaultBreakpointKmerSize,
                                                           echidna::caller::params::defaults::maxBreakpointKmerSize );

    // Truth data has 'event' right-aligned 10 bps due to not considering snp at 50664064.
    // Use these positions for convenience.
    const auto eventStart = 50664005;
    const auto eventEnd = 50664065;

    const auto breakpoint1 = std::make_shared< Breakpoint >(
        contig, eventStart + 20, true,
        referenceSequence->subseq( Region( contig, eventEnd + 20, eventEnd + 35 ) ).sequence() +
            BasePairSequence( "CAGAAAATGCCAGAAGGGAGGCC" ) );
    const auto breakpoint2 = std::make_shared< Breakpoint >(
        contig, eventEnd - 10, false,
        referenceSequence->subseq( Region( contig, eventStart - 42, eventStart - 10 ) ).sequence() );

    const auto nearSnp1 = std::make_shared< Variant >( referenceSequence, Region( contig, 50664064, 50664065 ), "G" );
    const auto nearSnp2 = std::make_shared< Variant >( referenceSequence, Region( contig, 50664094, 50664095 ), "G" );

    breakpoint2->addLocalVariants( {nearSnp1, nearSnp2} );

    const auto breakpoint3 = std::make_shared< Breakpoint >(
        contig, eventEnd - 10, false,
        referenceSequence->subseq( Region( contig, eventStart - 42, eventStart - 10 ) ).sequence() );

    auto startBpLocus = std::make_shared< BreakpointLocus >( contig, eventStart + 20, true );
    auto endBpLocus = std::make_shared< BreakpointLocus >( contig, eventEnd - 10, false );
    startBpLocus->add( breakpoint1 );
    endBpLocus->add( breakpoint2 );
    endBpLocus->add( breakpoint3 );

    ReadDataset readDataset( {""}, referenceSequence->region() );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint1, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint2, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint3, referenceSequence ) );
    const auto allReads = readDataset.getAllReads( 0 );

    const auto variants = breakpointVariantGenerator.getVariantCandidates( {startBpLocus, endBpLocus}, allReads );
    BOOST_REQUIRE_EQUAL( variants.size(), 2 );
    const std::vector< echidna::variant::varPtr_t > vecVariant( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( *vecVariant[0],
                       Variant( referenceSequence, Region( contig, eventStart - 10, eventEnd - 10 ), "" ) );
    BOOST_CHECK_EQUAL( *vecVariant[1], Variant( referenceSequence, Region( contig, eventStart, eventEnd ), "" ) );
}

BOOST_AUTO_TEST_CASE( testShouldFindDeletionAtChr1_165731209 )
{
    // 1000G known structural variant: DEL_pindel_1530
    const auto contig = "1";
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( contig, 165730859, 165731849 ),
        "ATGTCAAATTTCACATTTTAAAATGTTATAAAAATAACATTTTATGTCGGGTGCAGTGGCTCAGTTCTGTAAGCCCAGCACTCTGAGAGGCCAAGATGGGAAGATCATTT"
        "CAGTTCAGGAGTTCGAGATCAGCATGGGCAACACGGTGAAACCTCGTCTCTACCAAAACTACAAAAATCAGCCGGCCATGGTGGTGCACTCCTGTGGTCCCAGCTACTCA"
        "GGAGGCTGAGATGAGAGGATCATGGAAGCCCAGGAAGTCAAGACTGCAGTGAGCCATGATTGGGCCACTGCACTCCAACCTGGGTGACAGAGTGAAACCCTGTCTCAAAA"
        "AAATAAATAAATAAATAATAAACATTTTTATTAAAATGTTATAAGGATAACATTTTTATTAAAATGTTATAAGGATAACATTTTTATTAAAATGTTATAAGGATAACATT"
        "TTTATTAAAATGTTATAAGGATAACATTTTTATTAAAATGTTATAAGGATAACATTTTTATTAAAATGTTATAAGGATAACATTTTTATTAAAATGTTATAAGGATAACA"
        "TTTTTATTAAAATGTTATAAGGATAACATTTTTATTAAAATGTTATAAGGATAACATTTTTATTAAAATGTTATAAGGATAACATTTTTATTAAAATGTTATAAGGATAA"
        "CATTTTAAAAGACAAAATGGACAGAAAACACTTAATGTTTTGTTCTTTATCATTGTTATTAAACAAAGAAACAAAAATATTTCTTTCCCACCCCAGAAGGACTATAACAC"
        "AACTACTTGGCATTTAACCATGACTGCATAGAAAACAGGGCAATAAAAAGACAGAAGGCTAGAGAATAGAAACATAAGACAATTCCATAAACATTAGTCTGGTCTTCTAG"
        "AAGACCAGATTGAGACTAGATAAATTGTTTTTAGGGATAAATTAATCTAATGTTTTTCAAACTTCATGAGGTTGTGATGAAACCTTGGTGTGTGGATATTTATTCATCA"
        "G" );

    const auto eventStart = 165731209;
    const auto eventEnd = 165731489;
    const auto breakpoint1 = std::make_shared< Breakpoint >(
        contig, eventStart + 36, true,
        referenceSequence->subseq( Region( contig, eventEnd + 36, eventEnd + 72 ) ).sequence() );

    const auto breakpoint2 = std::make_shared< Breakpoint >(
        contig, eventEnd, false,
        referenceSequence->subseq( Region( contig, eventStart - 35, eventStart ) ).sequence() );

    auto startBpLocus = std::make_shared< BreakpointLocus >( contig, eventStart + 36, true );
    auto endBpLocus = std::make_shared< BreakpointLocus >( contig, eventEnd, false );

    ReadDataset readDataset( {""}, referenceSequence->region() );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint1, referenceSequence ) );
    readDataset.insertRead( "", readFromBreakpoint( *breakpoint2, referenceSequence ) );
    const auto allReads = readDataset.getAllReads( 0 );

    BreakpointVariantGenerator breakpointVariantGenerator( referenceSequence,
                                                           echidna::caller::params::defaults::defaultBreakpointKmerSize,
                                                           echidna::caller::params::defaults::maxBreakpointKmerSize );
    const auto variants = breakpointVariantGenerator.getVariantCandidates( {startBpLocus, endBpLocus}, allReads );

    BOOST_REQUIRE_EQUAL( variants.size(), 1 );
    const std::vector< echidna::variant::varPtr_t > vecVariant( variants.begin(), variants.end() );

    BOOST_CHECK_EQUAL( *vecVariant[0], Variant( referenceSequence, Region( contig, eventStart, eventEnd ), "" ) );
}