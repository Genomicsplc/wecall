#include "io/tabixFile.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "ioFixture.hpp"

BOOST_AUTO_TEST_CASE( testBasicTabixReading )
{
    std::string filepath =
        echidna::test::BUILD_PIPELINE_TEST_DATA + "/" + "vcf_comparison_test_data/output_original.vcf.gz";
    echidna::io::TabixFile tabixFile( filepath, filepath + ".tbi" );
    auto lines = tabixFile.fetch( echidna::caller::Region( "20", 0, 100000 ) );

    std::vector< std::string > expectedLines = {
        "20\t61098\t.\tC\tT\t233\tPASS\tFILTS=;PP=233.00;TC=33;TCR=19;TCF=14;VC=10;VCR=7;VCF=3;MQ=60.00;ABPV=0.05895;"
        "SBPV=0.8461\tGT:GL:GQ:NR:NV\t0/1:-27.06,0.30,-28.86:100.00:33:10"};

    BOOST_CHECK_EQUAL_COLLECTIONS( expectedLines.begin(), expectedLines.end(), lines.begin(), lines.end() );
}
