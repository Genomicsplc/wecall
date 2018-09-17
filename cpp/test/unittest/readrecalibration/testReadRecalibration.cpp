#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "readrecalibration/readRecalibration.hpp"
#include "io/read.hpp"
#include "alignment/cigar.hpp"

#include <memory>

using Read = echidna::io::Read;
using echidna::alignment::Cigar;
using echidna::alignment::cigarFlags;

// BOOST_AUTO_TEST_CASE( testFloorLowQualityScores )
//{
//    echidna::io::perSamReadRanges_t mapReadStartEnd;
//
//    char lowQual = 5;
//
//    echidna::io::readIntervalTree_t readContainer( 0, 100 );
//    readContainer.insert( std::make_shared< Read >( std::string( 3, 'T' ), std::string( 3, lowQual ), "", Cigar( "3M"
//    ),
//                                                    0, 0, 0, 1, 2, 2, 0 ) );
//
//    mapReadStartEnd["SomeSample"] = std::make_pair( readContainer.begin(), readContainer.end() );
//
//    echidna::corrector::floorLowQualityScores( mapReadStartEnd, 5, 2 );
//
//    for ( auto it = readContainer.begin(); it != readContainer.end(); ++it )
//    {
//        for ( const auto & qualityChar : it->getQualities() )
//        {
//            BOOST_CHECK_EQUAL( qualityChar, 2 );
//        }
//    }
//}
//
// BOOST_AUTO_TEST_CASE( testShouldNotFloorHighQualityScores )
//{
//    echidna::io::perSamReadRanges_t mapReadStartEnd;
//
//    char highQual = 5;
//
//    echidna::io::readIntervalTree_t readContainer( 0, 100 );
//    readContainer.insert( std::make_shared< Read >( std::string( 3, 'T' ), std::string( 3, highQual ), "",
//                                                    Cigar( "3M" ), 0, 0, 0, 1, 2, 2, 0 ) );
//
//    mapReadStartEnd["SomeSample"] = std::make_pair( readContainer.begin(), readContainer.end() );
//
//    echidna::corrector::floorLowQualityScores( mapReadStartEnd, 4, 2 );
//
//    for ( auto it = readContainer.begin(); it != readContainer.end(); ++it )
//    {
//        for ( const auto & qualityChar : it->getQualities() )
//        {
//            BOOST_CHECK_EQUAL( qualityChar, highQual );
//        }
//    }
//}
