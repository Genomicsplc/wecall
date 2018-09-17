#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "caller/region.hpp"
#include "caller/regionUtils.hpp"
#include "utils/flatten.hpp"
#include "ioTest/io/ioFixture.hpp"

using echidna::caller::Region;
using echidna::utils::Interval;
using echidna::caller::parseRegionString;

BOOST_FIXTURE_TEST_CASE( shouldExtractCorrectRegionsWhenNoneSuppliedOnCmdLine, echidna::test::FastaIndexFileFixture )
{
    auto actualRegions = echidna::utils::functional::flatten(
        echidna::caller::DataRegionsBuilder( {}, echidna::io::FastaIndex( indexFilename ) ).build() );

    // manual checks
    BOOST_CHECK_EQUAL( actualRegions.size(), 25 );
    BOOST_CHECK_EQUAL( actualRegions[0], Region( "1", 0, 249250621 ) );

    // automatic checks
    std::vector< Region > expectedRegions;
    for ( auto contigStr : fastaIndices[0]->standardContigs() )
    {
        expectedRegions.emplace_back( contigStr, fastaIndices[0]->contigs().at( contigStr ) );
    }
    BOOST_CHECK_EQUAL_COLLECTIONS( expectedRegions.begin(), expectedRegions.end(), actualRegions.begin(),
                                   actualRegions.end() );
}

BOOST_FIXTURE_TEST_CASE( shlouldExtractCorrectRegionsFromAList, echidna::test::FastaIndexFileFixture )
{
    auto actualRegions = echidna::utils::functional::flatten(
        echidna::caller::DataRegionsBuilder( {"1:2-3", "2:100-200"}, echidna::io::FastaIndex( indexFilename ) )
            .build() );
    std::vector< Region > expectedRegions = {Region( "1", 2, 3 ), Region( "2", 100, 200 )};

    BOOST_CHECK_EQUAL_COLLECTIONS( expectedRegions.begin(), expectedRegions.end(), actualRegions.begin(),
                                   actualRegions.end() );
}