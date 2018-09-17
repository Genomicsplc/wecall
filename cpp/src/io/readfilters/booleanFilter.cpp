// All content Copyright (C) 2018 Genomics plc
#include "io/readfilters/booleanFilter.hpp"
#include <sstream>
#include <iostream>

namespace echidna
{
namespace io
{
    //-----------------------------------------------------------------------------------------

    BooleanFilter::BooleanFilter( std::function< const bool(const Read &)> theFunction,
                                  const std::string & filterName,
                                  const bool negate )
        : m_function( theFunction ), m_filterName( filterName ), m_negate( negate )
    {
    }

    //-----------------------------------------------------------------------------------------

    bool BooleanFilter::passesFilter_impl( const io::Read & theRead )
    {
        if ( not this->m_negate )
        {
            return this->m_function( theRead );
        }
        else
        {
            return ( not this->m_function( theRead ) );
        }
    }

    //-----------------------------------------------------------------------------------------

    std::string BooleanFilter::toString() const
    {
        std::stringstream repr;
        repr << this->m_filterName << "()";
        return repr.str();
    }

    //-----------------------------------------------------------------------------------------
}
}
