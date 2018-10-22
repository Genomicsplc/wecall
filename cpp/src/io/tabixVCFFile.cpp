// All content Copyright (C) 2018 Genomics plc
#include <iostream>
#include <fstream>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/regex.hpp>
#include <utils/exceptions.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <variant/type/variant.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "io/tabixVCFFile.hpp"
#include "vcf/reader.hpp"
#include "vcf/filterDescription.hpp"

namespace wecall
{
namespace io
{

    TabixVCFFile::TabixVCFFile( std::string filename, std::string indexFilename )
        : m_tabixFile( filename, indexFilename )
    {
        this->readHeader( filename );
    }

    TabixVCFFile::~TabixVCFFile() {}

    void TabixVCFFile::readHeader( std::string filename )
    {
        ECHIDNA_LOG( DEBUG, "Reading vcf file header" );

        const auto headerLines = m_tabixFile.header();

        if ( headerLines.empty() )
        {
            throw utils::wecall_exception( "Malformed VCF: missing header." );
        }
        else if ( not boost::starts_with( headerLines.back(), "#CHROM" ) )
        {
            throw utils::wecall_exception( "Malformed VCF: did not find #CHROM at start of final header line." );
        }

        for ( const auto & line : headerLines )
        {
            if ( boost::starts_with( line, "##FILTER" ) )
            {
                ECHIDNA_LOG( DEBUG, "Found filter " << line );
                auto filterDesc = parseFilterHeaderLine( line );
                if ( containsFilterId( m_filterDescs, filterDesc.id ) )
                {
                    throw utils::wecall_exception( "The header contains repeated filter " + filterDesc.id + "." );
                }
                m_filterDescs.insert( filterDesc );
            }
            else if ( boost::starts_with( line, "##INFO" ) )
            {
                ECHIDNA_LOG( DEBUG, "VCF 'INFO' fields are not yet parsed. Skipping: " << line );
            }
            else if ( boost::starts_with( line, "##FORMAT" ) )
            {
                ECHIDNA_LOG( DEBUG, "VCF 'FORMAT' fields are not yet parsed. Skipping: " << line );
            }
            else if ( boost::starts_with( line, "##" ) )
            {
                auto metadata = parseMetaInfoLine( line );
                m_metaInformation.insert( metadata );
            }
            else if ( boost::starts_with( line, "#CHROM" ) )
            {
                // actual data on subsequent lines, the header has been parsed,
                // any additional header lines should be ignored as if they are comments.
                break;
            }
        }
    }

    std::vector< vcf::Record > TabixVCFFile::fetch( const caller::Region & region ) const
    {
        std::vector< vcf::Record > records;
        for ( const auto & line : m_tabixFile.fetch( region ) )
        {
            std::vector< std::string > cols;
            std::vector< std::string > alts;
            std::set< std::string > filterIds;

            boost::split( cols, line, boost::is_any_of( "\t" ) );

            if ( cols.size() < 8 )
            {
                throw utils::wecall_exception( "Record line should have at least 7 columns: " + line );
            }

            const std::string chrom = cols.at( 0 );
            const std::size_t one_based_pos = boost::lexical_cast< size_t >( cols.at( 1 ) );
            const auto stringIDs = cols.at( 2 );
            const std::string ref = cols.at( 3 );
            boost::split( alts, cols.at( 4 ), boost::is_any_of( "," ) );
            const double qual =
                cols.at( 5 ) == constants::vcfUnknownValue ? -1 : boost::lexical_cast< double >( cols.at( 5 ) );

            boost::split( filterIds, cols.at( 6 ), boost::is_any_of( ";" ) );

            const vcf::Info info = vcf::parseVCFInfo( cols.at( 7 ) );

            std::set< std::string > ids;
            if ( stringIDs != constants::vcfUnknownValue )
            {
                boost::split( ids, stringIDs, boost::is_any_of( ";" ) );
            }

            for ( const auto & filterId : filterIds )
            {
                if ( not containsFilterId( m_filterDescs, filterId ) and filterId != "PASS" )
                {
                    ECHIDNA_LOG( DEBUG, "Filter ID \"" + filterId + "\" was not supplied in the header." );
                }
            }

            vcf::SampleInfo emptySampleInfo;

            records.emplace_back( chrom, one_based_pos, ids, ref, alts, qual, filterIds, info, emptySampleInfo );
        }

        return records;
    }

    bool TabixVCFFile::containsFilterId( const std::set< vcf::FilterDesc > & filterDescs, const std::string & filterId )
    {
        // Try searching the set see if any filterDescs object contain the id
        auto it = std::find_if( filterDescs.begin(), filterDescs.end(), [&filterId]( const vcf::FilterDesc & desc )
                                {
                                    return desc.id == filterId;
                                } );
        return it != filterDescs.end();
    }

    std::pair< std::string, std::string > TabixVCFFile::parseMetaInfoLine( const std::string line )
    {
        boost::regex re( "^##(.+)=(.+)$" );
        boost::cmatch matches;
        if ( boost::regex_match( line.c_str(), matches, re ) )
        {
            const std::string key( matches[1].first, matches[1].second );
            const std::string value( matches[2].first, matches[2].second );
            return std::make_pair( key, value );
        }
        else
        {
            throw utils::wecall_exception( "Failed to match meta information header line: " + line );
        }
    }

    vcf::FilterDesc TabixVCFFile::parseFilterHeaderLine( const std::string & line )
    {
        boost::regex re( "^##FILTER=<ID=([^,]+),\\s*Description=\\\"(.*)\\\">\\s*$" );
        boost::cmatch matches;

        if ( boost::regex_match( line.c_str(), matches, re ) )
        {
            const std::string id( matches[1].first, matches[1].second );
            const std::string desc( matches[2].first, matches[2].second );
            return vcf::FilterDesc( id, desc );
        }
        else
        {
            throw utils::wecall_exception( "Failed to match filter header line: " + line );
        }
    }
}
}
