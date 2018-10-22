// All content Copyright (C) 2018 Genomics plc
#pragma once

namespace wecall
{
namespace utils
{
    namespace functional
    {

        template < typename T >
        T identity( T t )
        {
            return t;
        }

    }  // namespace functional
}  // namespace utils
}  // namespace wecall
