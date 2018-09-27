// All content Copyright (C) 2018 Genomics plc
#include "variant/type/variant.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <string>

using echidna::caller::Region;
using echidna::utils::ReferenceSequence;
using echidna::variant::Variant;
using echidna::variant::varPtr_t;

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( testSNP )
{
    std::string contig = "1";

    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( contig, 1000000, 1000001 ), "A" );
    varPtr_t theSnp = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "T" );

    BOOST_CHECK_EQUAL( theSnp->region(), Region( contig, 1000000, 1000001 ) );
    BOOST_CHECK_EQUAL( theSnp->refSequence().sequence(), "A" );
    BOOST_CHECK_EQUAL( theSnp->sequence(), "T" );
    BOOST_CHECK_EQUAL( theSnp->start(), 1000000 );
    BOOST_CHECK_EQUAL( theSnp->end(), 1000001 );

    auto expectedPrior = 1e-3 / 3.0;

    echidna::variant::setDefaultPriors( {theSnp} );
    BOOST_CHECK_CLOSE( theSnp->prior(), expectedPrior, 1e-5 );
}

BOOST_AUTO_TEST_CASE( testSnpLeftAlignDoesntChangePosition )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 9 ), "ATATATATA" );
    varPtr_t theSnp = std::make_shared< Variant >( referenceSequence, Region( "1", 8, 9 ), "T" );
    theSnp->getLeftAligned( 0 );
    BOOST_CHECK_EQUAL( theSnp->start(), 8 );
}
