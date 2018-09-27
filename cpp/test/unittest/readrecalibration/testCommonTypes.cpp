// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "readrecalibration/commonTypes.hpp"
#include "stats/functions.hpp"

BOOST_AUTO_TEST_CASE( testPhredToPCache )
{
    for ( auto i = 0; i < 100; ++i )
    {
        auto cacheNumber = echidna::corrector::phred_to_p( i );
        auto computedNumber = echidna::stats::fromPhredQ( i );
        BOOST_CHECK_CLOSE( cacheNumber, computedNumber, 1e-8 );
    }
}
