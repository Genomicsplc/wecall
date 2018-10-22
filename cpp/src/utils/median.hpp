// All content Copyright (C) 2018 Genomics plc
#pragma once

#include "utils/logging.hpp"

namespace wecall
{
namespace utils
{
    namespace functional
    {
        template < typename Type >
        Type median( const std::vector< Type > & values )
        {
            std::vector< Type > vecTypes = values;
            ECHIDNA_ASSERT( not vecTypes.empty(), "Median cannot be computed of empty list" );

            std::size_t middleIdx = vecTypes.size() / 2;

            auto target = vecTypes.begin() + middleIdx;
            std::nth_element( vecTypes.begin(), target, vecTypes.end() );

            if ( vecTypes.size() % 2 != 0 )
            {
                return *target;
            }
            else
            {
                auto targetNeighbour = std::max_element( vecTypes.begin(), target );
                return ( *target + *targetNeighbour ) / 2;
            }
        }
    }
}
}