// All content Copyright (C) 2018 Genomics plc
#include <algorithm>
#include "io/readUtils.hpp"

namespace wecall
{
namespace io
{
    namespace read
    {
        phred_t minBaseQualityInReadAroundInterval( const Read & read,
                                                    const utils::Interval & interval,
                                                    const int64_t padding )
        {
            const auto intervalInReadSpace = read.getIntervalInRead( interval );
            const auto paddedInterval = intervalInReadSpace.getPadded( padding );
            const auto intersectPadded = paddedInterval.getIntersect( utils::Interval( 0L, read.getLength() ) );

            const auto & qualities = read.getQualities();

            const auto min_element = std::min_element( qualities.cbegin() + intersectPadded.start(),
                                                       qualities.cbegin() + intersectPadded.end() );

            if ( min_element == qualities.cend() )
            {
                return 0;
            }
            else
            {
                return static_cast< phred_t >( *min_element );
            }
        }
    }
}
}
