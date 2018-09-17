// All content Copyright (C) 2018 Genomics plc
#ifndef BEST_SCORE_SELECTOR
#define BEST_SCORE_SELECTOR

#include <cstdint>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <algorithm>
#include <cassert>
#include <numeric>
#include "common.hpp"

namespace echidna
{
namespace utils
{
    std::size_t indexOfHighestValue( const std::vector< double > & values );

    template < typename T >
    std::vector< std::size_t > indiciesWithHighestValues( const std::vector< T > & values,
                                                          const std::size_t totalAllowed,
                                                          const double fractionOfTotalSumToConsider = 1.0,
                                                          const double maxHighestLowestRatio = 1.0e5 )
    {
        if ( values.empty() )
        {
            return {};
        }

        const T maxValue = *std::max_element( values.cbegin(), values.cend() );
        const T minValue = maxValue / maxHighestLowestRatio;

        std::vector< std::pair< T, std::size_t > > frequenciesToIndex;

        for ( auto it = values.cbegin(), end = values.cend(), begin = values.cbegin(); it != end; it++ )
        {
            if ( *it > minValue )
            {
                frequenciesToIndex.push_back( std::make_pair( *it, it - begin ) );
            }
        }

        std::sort( frequenciesToIndex.begin(), frequenciesToIndex.end() );

        const T totalCount = std::accumulate( values.cbegin(), values.cend(), T( 0 ) );
        const T countLimit = totalCount * fractionOfTotalSumToConsider;

        std::vector< std::size_t > bestIndicies;
        T cumulativeCount( 0 );
        auto frequenciesToIndexIt = frequenciesToIndex.rbegin();
        auto frequenciesToIndexEnd = frequenciesToIndex.rend();

        for ( std::size_t i = 0; i < totalAllowed; ++i )
        {
            if ( frequenciesToIndexIt == frequenciesToIndexEnd or cumulativeCount > countLimit )
            {
                break;
            }

            cumulativeCount += frequenciesToIndexIt->first;
            bestIndicies.push_back( frequenciesToIndexIt->second );
            frequenciesToIndexIt++;
        }

        return bestIndicies;
    }
}
}

#endif
