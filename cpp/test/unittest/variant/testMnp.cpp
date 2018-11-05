// All content Copyright (C) 2018 Genomics plc
#include "variant/type/variant.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <string>

using wecall::utils::ReferenceSequence;
using wecall::caller::Region;

BOOST_AUTO_TEST_CASE( mnp )
{
    using namespace wecall::variant;

    std::string contig = "1";
    const std::string removed( "ABCD" );
    const std::string added( "HELL" );

    const auto referenceSequence =
        std::make_shared< ReferenceSequence >( Region( contig, 1000000, 1000000 + removed.size() ), removed );
    varPtr_t theMNP = std::make_shared< Variant >( referenceSequence, referenceSequence->region(), added );

    BOOST_CHECK_EQUAL( theMNP->sequence(), added );
    BOOST_CHECK_EQUAL( theMNP->refSequence().sequence(), removed );

    BOOST_CHECK_EQUAL( theMNP->start(), 1000000 );
    BOOST_CHECK_EQUAL( theMNP->end(), 1000004 );

    auto nDiffs = 4;
    auto expectedPrior = 5e-5 * pow( 0.1, nDiffs - 1 ) * ( 1.0 - 0.1 );

    wecall::variant::setDefaultPriors( {theMNP} );
    BOOST_CHECK_CLOSE( theMNP->prior(), expectedPrior, 1e-5 );
}

BOOST_AUTO_TEST_CASE( testMnpLeftAlignDoesntChangePosition )
{
    using namespace wecall::variant;

    std::string contig = "1";
    const std::string seq = "TT";

    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( contig, 0, 10 ), "ATATATATAA" );

    varPtr_t theMnp = std::make_shared< Variant >( referenceSequence, Region( contig, 8, 10 ), seq );

    theMnp->getLeftAligned( 0 );

    BOOST_CHECK_EQUAL( theMnp->start(), 8 );
}
