// All content Copyright (C) 2018 Genomics plc
#pragma once

#include <vector>
#include <algorithm>
#include <numeric>

namespace echidna
{
namespace utils
{

    template < typename T, typename Iterator >
    T indexedProduct( const std::vector< T > & values, const Iterator index_begin, const Iterator index_end, T init )
    {
        std::vector< T > intermediate( static_cast< std::size_t >( std::distance( index_begin, index_end ) ) );
        const auto lookup = [&values]( const std::size_t hapIndex ) -> T
        {
            return values.at( hapIndex );
        };
        std::transform( index_begin, index_end, intermediate.begin(), lookup );
        return std::accumulate( intermediate.cbegin(), intermediate.cend(), init, std::multiplies< T >() );
    }
}
}