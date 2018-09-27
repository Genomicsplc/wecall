// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "io/readUtils.hpp"
#include "utils/interval.hpp"

using echidna::io::Read;
using echidna::caller::Region;
using echidna::alignment::Cigar;
using echidna::io::read::minBaseQualityInReadAroundInterval;
using echidna::utils::Interval;
using echidna::utils::ReferenceSequence;

std::string createQualityString( const std::vector< phred_t > & qualities )
{
    std::string qualityString( qualities.size(), 0 );
    for ( std::size_t index = 0; index < qualities.size(); ++index )
    {
        qualityString[index] = static_cast< char >( qualities[index] );
    }
    return qualityString;
}

BOOST_AUTO_TEST_CASE( shouldFindCorrectBasePairInFlatAlignedRead )
{
    auto refSequence = std::make_shared< ReferenceSequence >( Region( "1", -1, 10 ), std::string( 11, 'A' ) );
    auto testRead = std::make_shared< Read >( std::string( 5, 'A' ), createQualityString( {0, 1, 2, 3, 4} ), "",
                                              Cigar( "5M" ), 0, 0, 0, 0, 0, 0, 0, refSequence );
    const phred_t expectedResult = 2;
    BOOST_CHECK_EQUAL( expectedResult, minBaseQualityInReadAroundInterval( *testRead, Interval( 2, 3 ), 0 ) );
}

BOOST_AUTO_TEST_CASE( shouldFindCorrectBasePairInFlatAlignedReadWithPadding )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 10 ), std::string( 10, 'A' ) );
    auto testRead = std::make_shared< Read >( std::string( 5, 'A' ), createQualityString( {0, 1, 2, 3, 4} ), "",
                                              Cigar( "5M" ), 0, 0, 0, 0, 0, 0, 0, refSequence );
    const phred_t expectedResult = 1;
    BOOST_CHECK_EQUAL( expectedResult, minBaseQualityInReadAroundInterval( *testRead, Interval( 2, 3 ), 1 ) );
}

BOOST_AUTO_TEST_CASE( shouldGetValuesFromInsersionsAroundMatch )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 200 ), std::string( 200, 'A' ) );
    const Cigar cigar( "1M1I1M1I1M" );
    const auto startPos = 100L;

    auto testRead = std::make_shared< Read >( std::string( 5, 'A' ), createQualityString( {10, 50, 100, 51, 11} ), "",
                                              cigar, 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const phred_t expectedResult = 50;
    BOOST_CHECK_EQUAL( expectedResult,
                       minBaseQualityInReadAroundInterval( *testRead, Interval( startPos + 1, startPos + 2 ), 0 ) );
}
