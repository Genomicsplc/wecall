#include "variant/variantContainer.hpp"
#include "variant/type/variant.hpp"
#include "alignment/cigar.hpp"
#include "io/readDataSet.hpp"
#include "variant/haplotypeGenerator.hpp"
#include "unittest/vcf/VCFTestUtils.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <iostream>

using namespace echidna::alignment;
using echidna::utils::ReferenceSequence;
using echidna::caller::Region;
using echidna::variant::Variant;
using echidna::variant::varPtr_t;
using echidna::variant::VariantContainer;

echidna::io::readPtr_t make_read()
{
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( echidna::caller::Region( "1", 0, 100 ),
                                                                              std::string( 100, 'A' ) );
    return std::make_shared< echidna::io::Read >( "EDWARD", "ADRIAN", "", echidna::alignment::Cigar( "6M" ), 0, 0, 0, 0,
                                                  0, 0, 0, refSequence );
}

BOOST_AUTO_TEST_CASE( testVariantCountsComputesCorrectPercentCoverage )
{
    VariantContainer::VariantCounts variantCounts;
    variantCounts.m_totalReads = 0;
    variantCounts.m_totalVariantSupportingReads = 0;

    BOOST_CHECK_EQUAL( variantCounts.getPercentVariantCoverage(), 0 );

    variantCounts.m_totalReads = 3;
    variantCounts.m_totalVariantSupportingReads = 2;
    BOOST_CHECK_EQUAL( variantCounts.getPercentVariantCoverage(), 67 );
}

BOOST_AUTO_TEST_CASE( testRaisesIfNonSequentialVariantsAdded )
{

    VariantContainer variantContainer( 0, 0 );

    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    auto snp1 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "T" );
    auto snp2 = std::make_shared< Variant >( referenceSequence, Region( "1", 3, 4 ), "T" );
    auto read = make_read();

    BOOST_CHECK_THROW( variantContainer.addVariantsFromRead( read, {snp2, snp1}, {}, "sample" ),
                       echidna::utils::echidna_exception );
}

BOOST_AUTO_TEST_CASE( testShouldAddReadsToSingleVariant )
{

    VariantContainer variantContainer( 0, 0 );

    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    const Region varRegion = Region( "1", 1, 2 );
    const echidna::utils::BasePairSequence alt = "T";

    auto snp1 = std::make_shared< Variant >( referenceSequence, varRegion, alt );
    auto snp2 = std::make_shared< Variant >( referenceSequence, varRegion, alt );

    auto read1 = make_read();
    auto read2 = make_read();

    variantContainer.addVariantsFromRead( read1, {snp1}, {}, "sample" );
    variantContainer.addVariantsFromRead( read2, {snp2}, {}, "sample" );

    const auto variants = variantContainer.getVariants();

    BOOST_REQUIRE_EQUAL( 1, variants.size() );
    auto actualVar = *variants.begin();

    BOOST_CHECK_EQUAL( actualVar->region(), varRegion );
    BOOST_CHECK_EQUAL( actualVar->sequence(), alt );

    auto actualReads = actualVar->getReads();
    auto expectedReads = std::vector< echidna::io::readPtr_t >{read1, read2};
    BOOST_CHECK_EQUAL_COLLECTIONS( expectedReads.begin(), expectedReads.end(), actualReads.begin(), actualReads.end() );
}

BOOST_AUTO_TEST_CASE( testVariantContainerGetsCorrectSupport )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    std::vector< std::string > samples = {"Sample1, Sample2"};
    echidna::variant::VariantContainer variantContainer( 0, 0 );
    auto snp1 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "T" );
    auto snp2 = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "T" );
    variantContainer.addVariantsFromRead( make_read(), {snp1}, {}, samples.front() );
    variantContainer.addVariantsFromRead( make_read(), {snp2}, {}, samples.back() );

    const auto variants = variantContainer.getVariants();

    BOOST_CHECK_EQUAL( variants.size(), 1 );
    const std::vector< echidna::variant::varPtr_t > vecVariants( variants.begin(), variants.end() );
    BOOST_CHECK( checkVariantInVector( vecVariants,
                                       std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "T" ) ) );

    const auto count = variantContainer.totalReadsSupportingVariant( snp1 );
    BOOST_CHECK_EQUAL( count, 2 );

    auto newMatchingVariant = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "T" );
    const auto newCount = variantContainer.totalReadsSupportingVariant( newMatchingVariant );
    BOOST_CHECK_EQUAL( newCount, 2 );

    auto differentVariant = std::make_shared< Variant >( referenceSequence, Region( "2", 1, 2 ), "T" );
    BOOST_CHECK_THROW( variantContainer.totalReadsSupportingVariant( differentVariant ), std::out_of_range );

    std::shared_ptr< echidna::io::ReadDataset > readDataset =
        std::make_shared< echidna::io::ReadDataset >( samples, echidna::caller::Region( "1", 0, 10 ) );
    auto read = std::make_shared< echidna::io::Read >( std::string( 10, 'T' ), std::string( 10, 'Q' ), "test",
                                                       Cigar( "10M" ), 0, 0, 0, 0, 0, 0, 0, referenceSequence );
    readDataset->insertRead( samples.front(), read );
    readDataset->insertRead( samples.front(), read );
    variantContainer.computeCoverage( readDataset->region(), readDataset->getAllReads( 0 ) );

    const auto percent = variantContainer.maxReadPercentVariantCoverage( snp1 );
    BOOST_CHECK_EQUAL( percent, 100 );
}
//

BOOST_AUTO_TEST_CASE( testAddingCandidateVariantDoesntAddToSupport )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );
    echidna::variant::VariantContainer variantContainer( 0, 0 );
    auto snp = std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "T" );
    snp->disableFiltering();

    variantContainer.addCandidateVariant( snp, 1.0 );

    const auto variants = variantContainer.getVariants();

    BOOST_CHECK_EQUAL( variants.size(), 1 );
    const std::vector< echidna::variant::varPtr_t > vecVariants( variants.begin(), variants.end() );
    BOOST_CHECK( checkVariantInVector( vecVariants,
                                       std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "T" ) ) );

    const auto count = variantContainer.totalReadsSupportingVariant( snp );
    BOOST_CHECK_EQUAL( count, 0 );
}
