// All content Copyright (C) 2018 Genomics plc
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <numeric>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace echidna
{
namespace utils
{

    using matrix_t = boost::numeric::ublas::matrix< double >;
    using matrixRow_t = boost::numeric::ublas::matrix_row< const utils::matrix_t >;
    using matrixColumn_t = boost::numeric::ublas::matrix_column< const utils::matrix_t >;

    template < typename MatrixRow >
    double sumMatrixRowOverAllIndices( const MatrixRow & matrixRow )
    {
        return std::accumulate( matrixRow.begin(), matrixRow.end(), 0.0 );
    }

    template < typename MatrixRow, typename IndexSet >
    double sumMatrixRowOverIndexSubset( const MatrixRow & matrixRow, const IndexSet & indices )
    {
        auto theSum = 0.0;

        for ( const auto index : indices )
        {
            theSum += matrixRow( index );
        }
        return theSum;
    }

    void smoothLowOutliers( utils::matrix_t & matrix_t, const double maxDifference );

}  // namespace utils
}  // namespace echidna

#endif
