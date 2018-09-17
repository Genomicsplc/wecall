#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "variant/type/variant.hpp"
#include "caller/candidateVariantBank.hpp"
#include "unittest/vcf/VCFTestUtils.hpp"

using namespace echidna::variant;
using namespace echidna::caller;
using echidna::caller::Region;
using echidna::utils::ReferenceSequence;

BOOST_AUTO_TEST_CASE( shouldComputePriorFromAlleleFrequenceString )
{
    echidna::vcf::Info info;

    info.emplace_back( std::make_pair< std::string, std::vector< std::string > >(
        std::string( "AF" ), {std::string( "0.4" ), std::string( "0.5" )} ) );
    info.emplace_back( std::make_pair< std::string, std::vector< std::string > >(
        std::string( "DP" ), {std::string( "0.4" ), std::string( "0.1" )} ) );

    const auto expectedResults = {0.4, 0.5};
    const auto results = getPriorsFromInfo( info );

    BOOST_CHECK_EQUAL_COLLECTIONS( expectedResults.begin(), expectedResults.end(), results.begin(), results.end() );
}

//-------------------------------------------------------------------------------------------------
