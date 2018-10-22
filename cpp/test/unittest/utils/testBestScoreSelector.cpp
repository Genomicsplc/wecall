// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "utils/bestScoreSelector.hpp"

#include <vector>

BOOST_AUTO_TEST_CASE( testWithEmptySetOfScores )
{
    std::vector< double > scores = {};
    auto bestScores = wecall::utils::indiciesWithHighestValues( scores, 100 );
    BOOST_CHECK_EQUAL( bestScores.size(), 0 );
}

BOOST_AUTO_TEST_CASE( testGetsAllIndiciesIfTotalAllowedIsHigh )
{
    std::vector< double > scores = {1.0, 2.0};
    auto bestScores = wecall::utils::indiciesWithHighestValues( scores, 2 );
    BOOST_CHECK_EQUAL( bestScores.size(), 2 );

    BOOST_CHECK_EQUAL( bestScores[0], 1 );  // In reverse order as 2.0 > 1.0
    BOOST_CHECK_EQUAL( bestScores[1], 0 );
}

BOOST_AUTO_TEST_CASE( testGetsOnlyTheBestScore )
{
    std::vector< double > scores = {1.0, 2.0};
    auto bestScores = wecall::utils::indiciesWithHighestValues( scores, 1 );
    BOOST_CHECK_EQUAL( bestScores.size(), 1 );

    BOOST_CHECK_EQUAL( bestScores[0], 1 );
}

BOOST_AUTO_TEST_CASE( testShouldPickOutTheTwoHighestScores )
{
    std::vector< double > scores = {2.0, 1.0, 2.0};
    auto bestScores = wecall::utils::indiciesWithHighestValues( scores, 2 );
    BOOST_CHECK_EQUAL( bestScores.size(), 2 );

    BOOST_CHECK_EQUAL( scores[bestScores[0]], 2.0 );
    BOOST_CHECK_EQUAL( scores[bestScores[1]], 2.0 );
}

BOOST_AUTO_TEST_CASE( testShouldNotPickOutLowScores )
{
    std::vector< double > scores = {1.0e8, 1.0 + 1e-9, 1.0};
    const auto bestScores = wecall::utils::indiciesWithHighestValues( scores, scores.size(), 1.0, 1.0e8 );
    BOOST_CHECK_EQUAL( bestScores.size(), 2 );
    BOOST_CHECK_EQUAL( scores[bestScores[0]], 1e8 );
    BOOST_CHECK_EQUAL( scores[bestScores[1]], 1.0 + 1e-9 );
}
