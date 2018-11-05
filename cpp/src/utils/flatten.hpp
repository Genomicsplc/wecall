// All content Copyright (C) 2018 Genomics plc
#pragma once

namespace wecall
{
namespace utils
{
    namespace functional
    {

        template < template < typename... > class R = std::vector,
                   typename Outer,
                   typename Inner = typename Outer::value_type >
        R< typename Inner::value_type > flatten( const Outer & all )
        {
            R< typename Inner::value_type > flattened;

            for ( const auto & inner : all )
            {
                flattened.insert( std::end( flattened ), std::begin( inner ), std::end( inner ) );
            }

            return flattened;
        }

    }  // namespace functional
}  // namespace utils
}  // namespace wecall
