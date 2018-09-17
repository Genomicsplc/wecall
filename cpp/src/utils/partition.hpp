// All content Copyright (C) 2018 Genomics plc
#pragma once
#include <vector>

namespace echidna
{
namespace utils
{
    namespace functional
    {

        template < typename Input,
                   typename ComparisonFunction,
                   template < typename... > class R = std::vector,
                   typename Inner = Input >
        R< Inner > partition( const Input & all, const ComparisonFunction & comp )
        {
            R< Inner > partitioned;

            if ( not all.empty() )
            {
                partitioned.push_back( {*all.begin()} );

                for ( auto it = all.begin() + 1; it != all.end(); ++it )
                {
                    if ( comp( partitioned.back().back(), *it ) )
                    {
                        partitioned.back().push_back( *it );
                    }
                    else
                    {
                        partitioned.push_back( {*it} );
                    }
                }
            }

            return partitioned;
        }

    }  // namespace functional
}  // namespace utils
}  // namespace echidna
