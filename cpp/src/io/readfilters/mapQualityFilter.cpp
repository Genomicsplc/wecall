// All content Copyright (C) 2018 Genomics plc
#include "io/readfilters/mapQualityFilter.hpp"
#include <sstream>

namespace echidna
{
namespace io
{
    //-----------------------------------------------------------------------------------------

    bool MapQualFilter::passesFilter_impl( const io::Read & theRead )
    {
        return theRead.getMappingQuality() >= m_threshold;
    }

    //-----------------------------------------------------------------------------------------

    std::string MapQualFilter::toString() const
    {
        std::stringstream repr;
        repr << "MapQualFilter(Threshold = " << m_threshold << ")";
        return repr.str();
    }

    //-----------------------------------------------------------------------------------------
}
}
