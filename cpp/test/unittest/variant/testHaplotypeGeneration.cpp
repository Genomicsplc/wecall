// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <io/readDataSet.hpp>

#include "io/fastaFile.hpp"
#include "variant/type/variant.hpp"
#include "variant/haplotypeGenerator.hpp"

using echidna::caller::Region;
using echidna::utils::ReferenceSequence;
using echidna::utils::BasePairSequence;
using echidna::variant::AlignmentHaplotypeGenerator;
using echidna::variant::VariantCluster;
using echidna::variant::Variant;
using echidna::variant::varPtr_t;
using echidna::variant::variantSet_t;
using echidna::variant::varPtrComp;
using echidna::io::Read;
using echidna::alignment::Cigar;
using echidna::io::ReadDataset;
using echidna::utils::QualitySequence;
using echidna::io::RegionsReads;

std::shared_ptr< Read > makeFakeRead( echidna::utils::referenceSequencePtr_t reference,
                                      Region region,
                                      variantSet_t variants )
{
    const char defaultReadQual = 70;
    const auto readStartSeq = echidna::variant::Haplotype::buildHaplotypeSequence( reference, region, variants );
    const auto cigar = Cigar( std::to_string( region.size() ) + "M" );
    const auto mapq = 30;
    return std::make_shared< Read >( readStartSeq, QualitySequence( readStartSeq.size(), defaultReadQual ), "", cigar,
                                     0, region.start(), 0, mapq, 0, 0, 0, reference );
}

BOOST_AUTO_TEST_CASE( testHaplotypeGenerationWithLotsSNPsInTwoClustersAndUnsupported )
{

    const auto contig = "20";
    const auto referenceSequence = std::make_shared< ReferenceSequence >(
        Region( contig, 3037700, 3038250 ),
        "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        "TTATACTATGGTTCAAGACCAGCCTGGCCAACATGGGGAAACCCCATCTCTACTACAAAAATCAGCTGGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGG"
        "CCGAGGCGGGCGGATCACGAGGTCAAGAGATCGAGACCATCCCGGCTAAAACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAAATTAGCCGGGCGTAGTGGCGGGCG"
        "CCTGTAGTCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCCCGCCACTGCACTCCAGCCTGGGCGACAG"
        "AGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCAGCTGGGCATGATAGTATTACGTATAGGGGAGGGAATAAGAATATATATTANNNNNNNNN"
        "N" );

    const auto regionStart = 3037866;
    const auto regionEnd = 3038182;

    std::vector< varPtr_t > trueVariants = {
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 60, regionStart + 61 ), "A" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 61, regionStart + 62 ), "A" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 62, regionStart + 63 ), "A" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 63, regionStart + 64 ), "A" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 64, regionStart + 65 ), "T" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 65, regionStart + 66 ), "T" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 66, regionStart + 67 ), "A" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 120, regionStart + 121 ), "A" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 121, regionStart + 122 ), "A" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 122, regionStart + 123 ), "A" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 123, regionStart + 124 ), "C" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 124, regionStart + 125 ), "T" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 125, regionStart + 126 ), "A" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 126, regionStart + 127 ), "C" )};

    std::vector< varPtr_t > allVariants = {
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 5, regionStart + 6 ), "A" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 5, regionStart + 6 ), "G" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 5, regionStart + 6 ), "C" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 60, regionStart + 61 ), "T" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 60, regionStart + 61 ), "G" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 60, regionStart + 70 ), "" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 123, regionStart + 124 ), "G" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 124, regionStart + 125 ), "A" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 140, regionStart + 145 ), "" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 140, regionStart + 140 ),
                                     "TTTAT" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 200, regionStart + 201 ), "C" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 200, regionStart + 201 ), "T" ),
        std::make_shared< Variant >( referenceSequence, Region( contig, regionStart + 200, regionStart + 201 ), "G" ),
    };

    allVariants.insert( allVariants.end(), trueVariants.cbegin(), trueVariants.cend() );

    std::sort( allVariants.begin(), allVariants.end(), varPtrComp() );

    echidna::variant::setDefaultPriors( allVariants );

    VariantCluster variantCluster( allVariants, Region( contig, regionStart, regionEnd ) );

    ReadDataset readDataset( {"sample1", "sample2"}, referenceSequence->region() );
    readDataset.insertRead( "sample1",
                            makeFakeRead( referenceSequence, Region( contig, regionStart, regionStart + 100 ),
                                          variantSet_t( trueVariants.cbegin(), trueVariants.cbegin() + 7 ) ) );
    readDataset.insertRead( "sample1",
                            makeFakeRead( referenceSequence, Region( contig, regionStart + 40, regionStart + 140 ),
                                          variantSet_t( trueVariants.cbegin(), trueVariants.cend() ) ) );
    readDataset.insertRead( "sample1",
                            makeFakeRead( referenceSequence, Region( contig, regionStart + 110, regionStart + 210 ),
                                          variantSet_t( trueVariants.cbegin() + 7, trueVariants.cend() ) ) );

    int64_t maxHaplotypesAllowed( 2 );
    const auto reads = readDataset.getAllReads( 0 );
    AlignmentHaplotypeGenerator haplotypeGenerator( variantCluster.variants(), variantCluster.region(), reads,
                                                    referenceSequence, maxHaplotypesAllowed, 5 );
    const auto haps = haplotypeGenerator.generateHaplotypes();
    BOOST_REQUIRE_EQUAL( haps.size(), maxHaplotypesAllowed );
    BOOST_REQUIRE_EQUAL( haps[0].getVariants().size(), 14 );
    const auto actualVariants = haps[0].getVariants();
    const std::vector< varPtr_t > actualVecVars( actualVariants.cbegin(), actualVariants.cend() );

    for ( std::size_t i = 0; i < 14; ++i )
    {
        BOOST_CHECK_EQUAL( *trueVariants[i], *actualVecVars[i] );
    }

    BOOST_REQUIRE_EQUAL( haps[1].getVariants().size(), 0 );
}
