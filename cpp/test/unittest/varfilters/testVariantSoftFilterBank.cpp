#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/algorithm/string.hpp>

#include "varfilters/variantSoftFilterBank.hpp"
#include "utils/exceptions.hpp"
#include "vcf/filterDescription.hpp"

using VariantSoftFilterBank = echidna::varfilters::VariantSoftFilterBank;
using echidna_exception = echidna::utils::echidna_exception;

BOOST_AUTO_TEST_CASE( shouldConstructValidSoftFilterBank )
{
    std::vector< std::string > filters( {"AB", "SB"} );
    BOOST_CHECK_NO_THROW( VariantSoftFilterBank vsfb( filters, 0.01, 0.01, 0, 0, 0.0, 0.0, 0, 0 ) );
}

BOOST_AUTO_TEST_CASE( shouldRejectInvalidSoftFilter )
{
    std::vector< std::string > filters( {"ABBA"} );
    BOOST_CHECK_THROW( VariantSoftFilterBank vsfb( filters, 0.01, 0.01, 0, 0, 0.0, 0.0, 0, 0 ), echidna_exception );
}

BOOST_AUTO_TEST_CASE( shouldReturnFilterDescsInADeterministicOrder )
{
    // Note - Currently that deterministic order is lexicographically by filter ID
    std::vector< std::string > filters( {"SB", "AB"} );
    VariantSoftFilterBank vsfb( filters, 0.01, 0.01, 0, 0, 0.0, 0.0, 0, 0 );
    auto filterDescs = vsfb.getFilterDescs();
    BOOST_CHECK_EQUAL( filterDescs.size(), 2 );
    BOOST_CHECK_EQUAL( filterDescs[0].id, "AB" );
    BOOST_CHECK_EQUAL( filterDescs[1].id, "SB" );
}
