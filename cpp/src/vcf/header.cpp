// All content Copyright (C) 2018 Genomics plc
#include "vcf/header.hpp"

#include "utils/date.hpp"
#include "utils/write.hpp"

namespace echidna
{
namespace vcf
{
    Header::Header( const std::string & userSpecifiedFormat,
                    const caller::params::Application & applicationParams,
                    const std::string & ref,
                    const std::vector< std::string > & sampleNames,
                    const std::vector< std::string > & variantFieldIDs,
                    const std::vector< std::string > & genotypeFieldIDs,
                    const std::vector< FilterDesc > & filters,
                    const std::vector< caller::Region > & contigs )
        : m_userSpecifiedFormat( userSpecifiedFormat ),
          m_fileDate( utils::getCurrentDate() ),
          m_applicationParams( applicationParams ),
          m_ref( ref ),
          m_sampleNames( sampleNames ),
          m_infoFields( initInfoFields( variantFieldIDs ) ),
          m_formatFields( initFormatFields( genotypeFieldIDs ) ),
          m_filters( filters ),
          m_contigs( contigs )
    {
    }

    std::vector< Field > Header::initInfoFields( const std::vector< std::string > & variantFieldIDs )
    {
        std::vector< Field > fields;
        for ( const std::string & fieldID : variantFieldIDs )
        {
            if ( m_userSpecifiedFormat == constants::vcf41 )
            {
                fields.push_back( Field::infoFieldFromID( fieldID ) );
            }
            else
            {
                fields.push_back( Field::infoFieldFromID( fieldID, m_applicationParams.m_appName,
                                                          m_applicationParams.m_appVersion ) );
            }
        }
        return fields;
    }

    std::vector< Field > Header::initFormatFields( const std::vector< std::string > & genotypeFieldIDs )
    {
        std::vector< Field > fields;
        for ( const std::string & fieldID : genotypeFieldIDs )
        {
            if ( m_userSpecifiedFormat == constants::vcf41 )
            {
                fields.push_back( Field::formatFieldFromID( fieldID ) );
            }
            else
            {
                fields.push_back( Field::formatFieldFromID( fieldID, m_applicationParams.m_appName,
                                                            m_applicationParams.m_appVersion ) );
            }
        }
        return fields;
    }

    const std::map< std::string, std::string > Header::fileFormatMagicNumber = {{constants::vcf41, "VCFv4.1"},
                                                                                {constants::vcf42, "VCFv4.2"}};

    std::ostream & operator<<( std::ostream & out, const Header & header )
    {
        out << "##fileformat=" << Header::fileFormatMagicNumber.at( header.m_userSpecifiedFormat ) << std::endl;
        out << "##fileDate=" << header.m_fileDate << std::endl;
        out << "##source=" << header.m_applicationParams.m_appName << " v" << header.m_applicationParams.m_appVersion
            << std::endl;
        out << "##disclaimer=This software is in beta-testing. Results generated using the software are confidential "
               "and should only be used for research purposes in accordance with the legal agreement with Genomics plc."
            << std::endl;
        out << "##commit=" << header.m_applicationParams.m_appCommit << std::endl;
        out << "##buildTime=" << header.m_applicationParams.m_buildDate << std::endl;
        out << "##reference=" << header.m_ref << std::endl;
        out << "##options={" << header.m_applicationParams.m_appOptions << "}" << std::endl;

        for ( const auto & field : header.m_infoFields )
        {
            out << "##INFO=" << field << std::endl;
        }
        for ( const auto & filter : header.m_filters )
        {
            out << "##FILTER=" << filter << std::endl;
        }
        for ( const auto & field : header.m_formatFields )
        {
            out << "##FORMAT=" << field << std::endl;
        }
        for ( const auto & contig : header.m_contigs )
        {
            out << Header::contigKey() << "=<ID=" << contig.contig() << ",length=" << contig.end() << ">" << std::endl;
        }

        std::vector< std::string > columnHeaders = {
            Header::chromKey(), "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"};
        columnHeaders.insert( columnHeaders.end(), header.m_sampleNames.begin(), header.m_sampleNames.end() );
        out << boost::algorithm::join( columnHeaders, constants::vcfRecordColumnSeparator ) << std::endl;

        return out;
    }
}
}
