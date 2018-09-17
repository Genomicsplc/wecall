// All content Copyright (C) 2018 Genomics plc
#include "vcf/filterDescription.hpp"

namespace echidna
{
namespace vcf
{
    FilterDesc::FilterDesc( const std::string & id, const std::string & description )
        : id( id ), description( description )
    {
    }

    std::ostream & operator<<( std::ostream & out, const FilterDesc & filterDesc )
    {
        out << "<ID=" << filterDesc.id << ",Description=\"" << filterDesc.description << "\">";

        return out;
    }

    bool operator<( const FilterDesc & lhs, const FilterDesc & rhs ) { return ( lhs.id < rhs.id ); }
}
}
