#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "variant/type/variant.hpp"
#include "caller/callSet.hpp"

using echidna::caller::Call;
using echidna::variant::Variant;
using echidna::caller::Region;
using echidna::caller::CallComp;

// ---------------------------------------------------------------------------------------------------------------
// test: comparison of calls

BOOST_AUTO_TEST_CASE( shouldCompareVariantCalls )
{
    auto region1 = Region( "1", 1, 2 );
    auto region2 = Region( "1", 2, 3 );
    auto region3 = Region( "1", 3, 3 );
    auto region4 = Region( "2", 2, 3 );
    auto refSequenceChr1 =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'A' ) );
    auto refSequenceChr2 =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "2", 0, 5 ), std::string( 5, 'A' ) );

    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequenceChr1, region1, "C" );
    echidna::variant::varPtr_t variant2 = std::make_shared< Variant >( refSequenceChr1, region2, "T" );
    echidna::variant::varPtr_t variant3 = std::make_shared< Variant >( refSequenceChr1, region3, "CTA" );
    echidna::variant::varPtr_t variant4 = std::make_shared< Variant >( refSequenceChr2, region4, "C" );

    Call call1( variant1, variant1->interval(), 10, 1, {{Call::REF, Call::VAR}} );
    Call call2( variant2, variant2->interval(), 10, 1, {{Call::VAR, Call::VAR}} );
    Call call3( variant3, variant3->interval(), 10, 1, {{Call::REF, Call::VAR}} );
    Call call4( variant4, variant4->interval(), 10, 1, {{Call::VAR, Call::VAR}} );

    CallComp comp;

    BOOST_CHECK( comp( call1, call2 ) );
    BOOST_CHECK( comp( call2, call3 ) );
    BOOST_CHECK( comp( call3, call4 ) );
}

BOOST_AUTO_TEST_CASE( shouldCompareVariantAndRefCall )
{
    auto regionVar = Region( "1", 2, 3 );
    auto regionRef = Region( "1", 3, 4 );
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'A' ) );
    echidna::variant::varPtr_t variant1 = std::make_shared< Variant >( refSequence, regionVar, "C" );

    Call callVar( variant1, regionVar.interval(), 10, 1, {{Call::REF, Call::VAR}} );
    Call callRef( nullptr, regionRef.interval(), 10, 1, {{Call::REF, Call::REF}} );

    CallComp comp;

    BOOST_CHECK( comp( callVar, callRef ) );
}

BOOST_AUTO_TEST_CASE( shouldCompareRefCalls )
{
    // there should never be any overlapping ref calls
    auto region1 = Region( "1", 2, 5 );
    auto region2 = Region( "1", 5, 10 );
    auto region3 = Region( "1", 5, 20 );

    Call call1( nullptr, region1.interval(), 10, 1, {{Call::REF, Call::REF}} );
    Call call2( nullptr, region2.interval(), 10, 1, {{Call::REF, Call::REF}} );
    Call call3( nullptr, region3.interval(), 10, 1, {{Call::REF, Call::REF}} );

    CallComp comp;

    BOOST_CHECK( comp( call1, call2 ) );
    BOOST_CHECK( comp( call2, call3 ) );
}
