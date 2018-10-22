// All content Copyright (C) 2018 Genomics plc
#include "io/readDataSet.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <string>

#include "ioFixture.hpp"
#include "caller/params.hpp"
#include "caller/params.hpp"
#include "io/readDataSet.hpp"

using wecall::io::ReadDataset;

BOOST_AUTO_TEST_CASE( testReadDataSetInitialisation )
{
    std::vector< std::string > samples = {"NA12878", "NA12891"};
    ReadDataset readDataSet( samples, wecall::caller::Region( "1", 0, 0 ) );
    BOOST_CHECK( readDataSet.isEmpty() );
    std::vector< std::string > setSamples = readDataSet.getSampleNames();
    BOOST_CHECK_EQUAL_COLLECTIONS( setSamples.begin(), setSamples.end(), samples.begin(), samples.end() );
    auto readRange_NA12878 = readDataSet.getAllReads( 0 ).at( "NA12878" );
    BOOST_CHECK_EQUAL( std::distance( readRange_NA12878.begin(), readRange_NA12878.end() ), 0 );  // has 0 reads

    auto readRange_NA12891 = readDataSet.getAllReads( 0 ).at( "NA12891" );
    BOOST_CHECK_EQUAL( std::distance( readRange_NA12891.begin(), readRange_NA12891.end() ), 0 );  // has 0 reads
}
