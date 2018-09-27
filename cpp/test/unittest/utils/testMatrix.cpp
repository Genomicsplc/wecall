// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/matrix.hpp"
#include <set>
#include <cassert>
#include <algorithm>

echidna::utils::matrix_t getMatrixFromVecOfVecs( const std::vector< std::vector< double > > & values )
{
    assert( values.size() > 0 );

    echidna::utils::matrix_t matrix_t( values.size(), values[0].size() );

    for ( std::size_t rowIndex = 0; rowIndex < matrix_t.size1(); ++rowIndex )
    {
        const auto & row = values[rowIndex];
        assert( row.size() == matrix_t.size2() );

        std::copy( row.begin(), row.end(), matrix_t.data().begin() + matrix_t.size2() * rowIndex );
    }
    return matrix_t;
}

BOOST_AUTO_TEST_CASE( testSumOverMatrixRowIndexSet )
{
    std::vector< int > values = {1, 2, 3, 4, 5, 6, 7};
    echidna::utils::matrix_t matrix_t( 1, 7 );
    std::copy( values.begin(), values.end(), matrix_t.data().begin() );

    echidna::utils::matrixRow_t matrixRow_t( matrix_t, 0 );

    std::set< std::size_t > indicies = {0, 2, 4};
    BOOST_CHECK_EQUAL( echidna::utils::sumMatrixRowOverIndexSubset( matrixRow_t, indicies ), 1 + 3 + 5 );
    BOOST_CHECK_EQUAL( echidna::utils::sumMatrixRowOverAllIndices( matrixRow_t ), 1 + 2 + 3 + 4 + 5 + 6 + 7 );
}

BOOST_AUTO_TEST_CASE( testAdjustmentToMedian )
{
    std::vector< std::vector< double > > values = {
        {13.0}, {1.0}, {1.0}, {1.0e-6}, {1.0}, {1.0}, {1.0}, {15.0},
    };

    echidna::utils::matrix_t matrix = getMatrixFromVecOfVecs( values );

    BOOST_CHECK_CLOSE( *std::min_element( matrix.data().begin(), matrix.data().end() ), 1.0e-6, 1.0 );

    echidna::utils::smoothLowOutliers( matrix, 1.0e-4 );

    BOOST_CHECK_CLOSE( *std::min_element( matrix.data().begin(), matrix.data().end() ), 1.0e-4, 1.0 );
}