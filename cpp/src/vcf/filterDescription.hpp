// All content Copyright (C) 2018 Genomics plc
#ifndef VCF_FILTERDESC_HPP
#define VCF_FILTERDESC_HPP

#include "common.hpp"

#include <iostream>

namespace wecall
{
namespace vcf
{
    struct FilterDesc
    {
        /// Basic constructor from constituent data
        ///
        /// @param id Filter ID.
        /// @param description Filter description.
        FilterDesc( const std::string & id, const std::string & description );

        /// Writes the filter to the output stream in the form of a VCF FILTER header line.
        ///
        /// @param out Output stream
        /// @param filterDesc Filter to be output
        /// @return Output stream
        friend std::ostream & operator<<( std::ostream & out, const FilterDesc & filterDesc );
        friend bool operator<( const FilterDesc & lhs, const FilterDesc & rhs );

        std::string id;
        std::string description;
    };
}
}

#endif
