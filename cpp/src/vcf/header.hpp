// All content Copyright (C) 2018 Genomics plc
#ifndef VCF_HEADER_HPP
#define VCF_HEADER_HPP

#include "common.hpp"
#include "vcf/field.hpp"
#include "vcf/filterDescription.hpp"
#include "caller/params.hpp"
#include "caller/region.hpp"

#include <vector>
#include <string>
#include <iostream>

namespace wecall
{
namespace vcf
{
    /// An in-memory representation of a VCF header.
    class Header
    {
    public:
        /// Basic constructor from constituent data.
        ///
        /// @param source Name of the application generating the VCF output.
        /// @param version Version of the application generating the VCF output.
        /// @param commit The label associated with the application's current commit.
        /// @param ref Name of the reference file name.
        /// @param options Options selected for the application generating the VCF output.
        /// @param sampleNames Names of the samples from which variants were called.
        /// @param variantFieldIDs IDs of the INFO field definitions included in the header.
        /// @param genotypeFieldIDs IDs of the FORMAT field definitions included in the header.
        /// @param filters list of FILTERs whose definitions are to be included in the header.
        Header( const std::string & userSpecifiedFormat,
                const caller::params::Application & applicationParams,
                const std::string & ref,
                const std::vector< std::string > & sampleNames,
                const std::vector< std::string > & variantFieldIDs,
                const std::vector< std::string > & genotypeFieldIDs,
                const std::vector< FilterDesc > & filters,
                const std::vector< caller::Region > & contigs );

        /// Writes the VCF header to the output stream.
        ///
        /// @param out Output stream
        /// @param header Header to be output
        /// @return Output stream
        friend std::ostream & operator<<( std::ostream & out, const Header & header );

        // Slightly annoyingly contig has to be lower case according to vcf-validator.
        static std::string contigKey() { return "##contig"; }
        static std::string chromKey() { return "#CHROM"; }

    private:
        const static std::map< std::string, std::string > fileFormatMagicNumber;
        /// Initialise the header's INFO field definitions.
        ///
        /// @param variantFieldIDs Field IDs
        /// @return Field definitions
        std::vector< Field > initInfoFields( const std::vector< std::string > & variantFieldIDs );

        /// Initialise the header's FORMAT field definitions.
        ///
        /// @param genotypeFieldIDs Field IDs
        /// @return Field definitions
        std::vector< Field > initFormatFields( const std::vector< std::string > & genotypeFieldIDs );

        const std::string m_userSpecifiedFormat;
        const std::string m_fileDate;
        const caller::params::Application m_applicationParams;
        const std::string m_ref;
        const std::vector< std::string > m_sampleNames;
        const std::vector< Field > m_infoFields;
        const std::vector< Field > m_formatFields;
        const std::vector< FilterDesc > m_filters;
        const std::vector< caller::Region > m_contigs;
    };
}
}

#endif
