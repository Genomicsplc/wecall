// All content Copyright (C) 2018 Genomics plc
#include "utils/matrix.hpp"
#include "utils/median.hpp"

namespace wecall
{
namespace utils
{
    void smoothLowOutliers( utils::matrix_t & matrix_t, const double maxDifference )
    {
        if ( matrix_t.size1() == 0 or matrix_t.size2() == 0 )
        {
            return;
        }

        const auto nReads = matrix_t.size1();

        std::vector< double > maxValuesPerRead;

        for ( std::size_t readIndex = 0; readIndex < nReads; ++readIndex )
        {
            const utils::matrixRow_t row_t( matrix_t, readIndex );
            const auto max = *std::max_element( row_t.begin(), row_t.end() );
            maxValuesPerRead.push_back( max );
        }

        const auto median = utils::functional::median( maxValuesPerRead );

        const double minValue = median * maxDifference;

        for ( auto & val : matrix_t.data() )
        {
            val = std::max( val, minValue );
        }
    }
}
}