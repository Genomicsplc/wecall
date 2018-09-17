// All content Copyright (C) 2018 Genomics plc
#include "io/readfilters/baseQualityFilter.hpp"
#include <sstream>

namespace echidna
{
namespace io
{
    //-----------------------------------------------------------------------------------------

    bool BaseQualFilter::passesFilter_impl( const io::Read & theRead )
    {
        const auto & quals = theRead.getQualities();
        int count = 0;

        for ( int theQual : quals )
        {
            if ( theQual > this->m_qualThreshold )
            {
                ++count;

                if ( count > this->m_minBases )
                {
                    return true;
                }
            }
        }

        return false;
    }

    //-----------------------------------------------------------------------------------------

    std::string BaseQualFilter::toString() const
    {
        std::stringstream repr;
        repr << "BaseQualFilter(Threshold = " << m_qualThreshold << ", MinBases = " << m_minBases << ")";
        return repr.str();
    }

    //-----------------------------------------------------------------------------------------
}
}
