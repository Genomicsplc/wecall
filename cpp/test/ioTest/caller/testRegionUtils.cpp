// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "caller/region.hpp"
#include "caller/regionUtils.hpp"
#include "utils/flatten.hpp"
#include "ioTest/io/ioFixture.hpp"

using wecall::caller::Region;
using wecall::utils::Interval;
using wecall::caller::parseRegionString;

BOOST_FIXTURE_TEST_CASE( shouldExtractCorrectRegionsWhenNoneSuppliedOnCmdLine, wecall::test::FastaIndexFileFixture )
{
    auto actualRegions = wecall::utils::functional::flatten(
        wecall::caller::DataRegionsBuilder( {}, wecall::io::FastaIndex( indexFilename ) ).build() );

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

BOOST_FIXTURE_TEST_CASE( shlouldExtractCorrectRegionsFromAList, wecall::test::FastaIndexFileFixture )
{
    auto actualRegions = wecall::utils::functional::flatten(
        wecall::caller::DataRegionsBuilder( {"1:2-3", "2:100-200"}, wecall::io::FastaIndex( indexFilename ) )
            .build() );
    std::vector< Region > expectedRegions = {Region( "1", 2, 3 ), Region( "2", 100, 200 )};

    BOOST_CHECK_EQUAL_COLLECTIONS( expectedRegions.begin(), expectedRegions.end(), actualRegions.begin(),
                                   actualRegions.end() );
}