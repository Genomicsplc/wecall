// All content Copyright (C) 2018 Genomics plc
#include <cstdint>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <algorithm>
#include <cassert>
#include "common.hpp"
#include "utils/logging.hpp"
#include "utils/bestScoreSelector.hpp"

namespace wecall
{
namespace utils
{
    std::size_t indexOfHighestValue( const std::vector< double > & values )
    {
        return long_to_sizet( std::distance( values.cbegin(), std::max_element( values.cbegin(), values.cend() ) ) );
    }
}
}