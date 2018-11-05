// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "utils/NeedlemanWunsch.hpp"

using namespace wecall::variant;
using namespace wecall::caller;
using wecall::utils::NWVariant;
using wecall::utils::ReferenceSequence;

BOOST_AUTO_TEST_CASE( shouldProduceCorrectVariantForSNPAtStart )
{
    auto region = Region( "1", 10, 15 );
    auto ref_bases = "ATCAT";
    auto subRegion = Region( "1", 10, 11 );
    auto alt_bases = "G";

    auto referenceSequence = std::make_shared< ReferenceSequence >( region, ref_bases );
    auto originalVariant = Variant( referenceSequence, subRegion, alt_bases );

    int64_t base = region.start();
    auto nWVariant = NWVariant( subRegion.start() - base, subRegion.end() - base, alt_bases );

    BOOST_CHECK_EQUAL( originalVariant, *nWVariant.getVariant( region.contig(), base, referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( shouldProduceCorrectVariantForSNPInMiddle )
{
    auto region = Region( "1", 10, 15 );
    auto ref_bases = "ATCAT";
    auto subRegion = Region( "1", 11, 12 );
    auto alt_bases = "G";

    auto referenceSequence = std::make_shared< ReferenceSequence >( region, ref_bases );
    auto originalVariant = Variant( referenceSequence, subRegion, alt_bases );

    int64_t base = region.start();
    auto nWVariant = NWVariant( subRegion.start() - base, subRegion.end() - base, alt_bases );

    BOOST_CHECK_EQUAL( originalVariant, *nWVariant.getVariant( region.contig(), base, referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( shouldProduceCorrectVariantForSNPAtEnd )
{
    auto region = Region( "1", 10, 15 );
    auto ref_bases = "ATCAT";
    auto subRegion = Region( "1", 14, 15 );
    auto alt_bases = "G";

    auto referenceSequence = std::make_shared< ReferenceSequence >( region, ref_bases );
    auto originalVariant = Variant( referenceSequence, subRegion, alt_bases );

    int64_t base = region.start();
    auto nWVariant = NWVariant( subRegion.start() - base, subRegion.end() - base, alt_bases );

    BOOST_CHECK_EQUAL( originalVariant, *nWVariant.getVariant( region.contig(), base, referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( shouldProduceCorrectVariantForSingleBaseInsertionAtStart )
{
    auto region = Region( "1", 10, 15 );
    auto ref_bases = "ATCAT";
    auto subRegion = Region( "1", 10, 10 );
    auto alt_bases = "G";

    auto referenceSequence = std::make_shared< ReferenceSequence >( region, ref_bases );
    auto originalVariant = Variant( referenceSequence, subRegion, alt_bases );

    int64_t base = region.start();
    auto nWVariant = NWVariant( subRegion.start() - base, subRegion.end() - base, alt_bases );

    BOOST_CHECK_EQUAL( originalVariant, *nWVariant.getVariant( region.contig(), base, referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( shouldProduceCorrectVariantForSingleBaseInsertionInMiddle )
{
    auto region = Region( "1", 10, 15 );
    auto ref_bases = "ATCAT";
    auto subRegion = Region( "1", 12, 12 );
    auto alt_bases = "G";

    auto referenceSequence = std::make_shared< ReferenceSequence >( region, ref_bases );
    auto originalVariant = Variant( referenceSequence, subRegion, alt_bases );

    int64_t base = region.start();
    auto nWVariant = NWVariant( subRegion.start() - base, subRegion.end() - base, alt_bases );

    BOOST_CHECK_EQUAL( originalVariant, *nWVariant.getVariant( region.contig(), base, referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( shouldProduceCorrectVariantForSingleBaseInsertionAtEnd )
{
    auto region = Region( "1", 10, 15 );
    auto ref_bases = "ATCAT";
    auto subRegion = Region( "1", 15, 15 );
    auto alt_bases = "G";

    auto referenceSequence = std::make_shared< ReferenceSequence >( region, ref_bases );
    auto originalVariant = Variant( referenceSequence, subRegion, alt_bases );

    int64_t base = region.start();
    auto nWVariant = NWVariant( subRegion.start() - base, subRegion.end() - base, alt_bases );

    BOOST_CHECK_EQUAL( originalVariant, *nWVariant.getVariant( region.contig(), base, referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( shouldProduceCorrectVariantForSingleBaseDeletionAtStart )
{
    auto region = Region( "1", 10, 15 );
    auto ref_bases = "ATCAT";
    auto subRegion = Region( "1", 10, 11 );
    auto alt_bases = "";

    auto referenceSequence = std::make_shared< ReferenceSequence >( region, ref_bases );
    auto originalVariant = Variant( referenceSequence, subRegion, alt_bases );

    int64_t base = region.start();
    auto nWVariant = NWVariant( subRegion.start() - base, subRegion.end() - base, alt_bases );

    BOOST_CHECK_EQUAL( originalVariant, *nWVariant.getVariant( region.contig(), base, referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( shouldProduceCorrectVariantForSingleBaseDeletionInMiddle )
{
    auto region = Region( "1", 10, 15 );
    auto ref_bases = "ATCAT";
    auto subRegion = Region( "1", 12, 13 );
    auto alt_bases = "";

    auto referenceSequence = std::make_shared< ReferenceSequence >( region, ref_bases );
    auto originalVariant = Variant( referenceSequence, subRegion, alt_bases );

    int64_t base = region.start();
    auto nWVariant = NWVariant( subRegion.start() - base, subRegion.end() - base, alt_bases );

    BOOST_CHECK_EQUAL( originalVariant, *nWVariant.getVariant( region.contig(), base, referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( shouldProduceCorrectVariantForSingleBaseDeletionAtEnd )
{
    auto region = Region( "1", 10, 15 );
    auto ref_bases = "ATCAT";
    auto subRegion = Region( "1", 14, 15 );
    auto alt_bases = "";

    auto referenceSequence = std::make_shared< ReferenceSequence >( region, ref_bases );
    auto originalVariant = Variant( referenceSequence, subRegion, alt_bases );

    int64_t base = region.start();
    auto nWVariant = NWVariant( subRegion.start() - base, subRegion.end() - base, alt_bases );

    BOOST_CHECK_EQUAL( originalVariant, *nWVariant.getVariant( region.contig(), base, referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( shouldProduceCorrectVariantForMultiBaseInsertionAtStart )
{
    auto region = Region( "1", 10, 15 );
    auto ref_bases = "ATCAT";
    auto subRegion = Region( "1", 10, 10 );
    auto alt_bases = "GGGG";

    auto referenceSequence = std::make_shared< ReferenceSequence >( region, ref_bases );
    auto originalVariant = Variant( referenceSequence, subRegion, alt_bases );

    int64_t base = region.start();
    auto nWVariant = NWVariant( subRegion.start() - base, subRegion.end() - base, alt_bases );

    BOOST_CHECK_EQUAL( originalVariant, *nWVariant.getVariant( region.contig(), base, referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( shouldProduceCorrectVariantForMultiBaseInsertionInMiddle )
{
    auto region = Region( "1", 10, 15 );
    auto ref_bases = "ATCAT";
    auto subRegion = Region( "1", 12, 12 );
    auto alt_bases = "GGGG";

    auto referenceSequence = std::make_shared< ReferenceSequence >( region, ref_bases );
    auto originalVariant = Variant( referenceSequence, subRegion, alt_bases );

    int64_t base = region.start();
    auto nWVariant = NWVariant( subRegion.start() - base, subRegion.end() - base, alt_bases );

    BOOST_CHECK_EQUAL( originalVariant, *nWVariant.getVariant( region.contig(), base, referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( shouldProduceCorrectVariantForMultiBaseInsertionAtEnd )
{
    auto region = Region( "1", 10, 15 );
    auto ref_bases = "ATCAT";
    auto subRegion = Region( "1", 15, 15 );
    auto alt_bases = "GGGG";

    auto referenceSequence = std::make_shared< ReferenceSequence >( region, ref_bases );
    auto originalVariant = Variant( referenceSequence, subRegion, alt_bases );

    int64_t base = region.start();
    auto nWVariant = NWVariant( subRegion.start() - base, subRegion.end() - base, alt_bases );

    BOOST_CHECK_EQUAL( originalVariant, *nWVariant.getVariant( region.contig(), base, referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( shouldProduceCorrectVariantForMultiBaseDeletionAtStart )
{
    auto region = Region( "1", 10, 15 );
    auto ref_bases = "ATCAT";
    auto subRegion = Region( "1", 10, 12 );
    auto alt_bases = "";

    auto referenceSequence = std::make_shared< ReferenceSequence >( region, ref_bases );
    auto originalVariant = Variant( referenceSequence, subRegion, alt_bases );

    int64_t base = region.start();
    auto nWVariant = NWVariant( subRegion.start() - base, subRegion.end() - base, alt_bases );

    BOOST_CHECK_EQUAL( originalVariant, *nWVariant.getVariant( region.contig(), base, referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( shouldProduceCorrectVariantForMultiBaseDeletionInMiddle )
{
    auto region = Region( "1", 10, 15 );
    auto ref_bases = "ATCAT";
    auto subRegion = Region( "1", 12, 14 );
    auto alt_bases = "";

    auto referenceSequence = std::make_shared< ReferenceSequence >( region, ref_bases );
    auto originalVariant = Variant( referenceSequence, subRegion, alt_bases );

    int64_t base = region.start();
    auto nWVariant = NWVariant( subRegion.start() - base, subRegion.end() - base, alt_bases );

    BOOST_CHECK_EQUAL( originalVariant, *nWVariant.getVariant( region.contig(), base, referenceSequence ) );
}

BOOST_AUTO_TEST_CASE( shouldProduceCorrectVariantForMultiBaseDeletionAtEnd )
{
    auto region = Region( "1", 10, 15 );
    auto ref_bases = "ATCAT";
    auto subRegion = Region( "1", 13, 15 );
    auto alt_bases = "";

    auto referenceSequence = std::make_shared< ReferenceSequence >( region, ref_bases );
    auto originalVariant = Variant( referenceSequence, subRegion, alt_bases );

    int64_t base = region.start();
    auto nWVariant = NWVariant( subRegion.start() - base, subRegion.end() - base, alt_bases );

    BOOST_CHECK_EQUAL( originalVariant, *nWVariant.getVariant( region.contig(), base, referenceSequence ) );
}
