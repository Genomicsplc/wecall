// All content Copyright (C) 2018 Genomics plc
#include <fstream>
#include <iostream>
#include <vector>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem/path.hpp>
#include "io/bedFile.hpp"
#include "caller/region.hpp"
#include "caller/regionUtils.hpp"
#include "utils/exceptions.hpp"
#include "utils/logging.hpp"
#include "utils/timer.hpp"

namespace wecall
{
namespace io
{
    BedFile::BedFile( std::string fileName )
        : m_filename( fileName ), m_timer( std::make_shared< utils::Timer >( "IO", utils::fileMetaData( fileName ) ) )
    {
        WECALL_ASSERT( boost::filesystem::exists( m_filename ), "BED file " + m_filename + " does not exist" );
    }

    caller::regions_t BedFile::getRegions() const
    {
        caller::regions_t regions;
        if ( boost::filesystem::path( m_filename ).extension().string() == ".bed" )
        {
            std::ifstream inStream( m_filename, std::ios_base::in );
            regions = this->readRegionsFromStream( inStream );
            inStream.close();
        }
        else if ( boost::filesystem::path( m_filename ).extension().string() == ".gz" )
        {
            std::ifstream file( m_filename, std::ios_base::in );

            boost::iostreams::filtering_streambuf< boost::iostreams::input > inbuf;
            inbuf.push( boost::iostreams::gzip_decompressor() );
            inbuf.push( file );

            std::istream inStream( &inbuf );
            regions = this->readRegionsFromStream( inStream );

            file.close();
        }

        return regions;
    }

    caller::regions_t BedFile::readRegionsFromStream( std::istream & inStream ) const
    {
        utils::ScopedTimerTrigger scopedTimerTrigger( m_timer );
        caller::regions_t regions;

        std::string line;
        size_t lineNumber = 0;
        bool insideHeader = true;

        while ( std::getline( inStream, line ) )
        {
            ++lineNumber;

            auto columns = parseRegionLine( line, insideHeader, lineNumber );
            if ( not columns.empty() )
            {
                insideHeader = false;

                regions.emplace_back( boost::lexical_cast< std::string >( columns[0] ),
                                      boost::lexical_cast< int64_t >( columns[1] ),
                                      boost::lexical_cast< int64_t >( columns[2] ) );
            }
        }

        return regions;
    }

    std::vector< std::string > BedFile::parseRegionLine( std::string line, bool insideHeader, size_t lineNumber ) const
    {
        boost::regex bedLineStart( "^" + caller::chromMatch + "\\t" + caller::posMatch + "\\t" + caller::posMatch +
                                   ".*" );
        boost::cmatch what;
        if ( boost::regex_match( line.c_str(), what, bedLineStart ) )
        {
            return {what[1], what[2], what[3]};
        }
        else if ( line.empty() or insideHeader )
        {
            return {};
        }
        else
        {
            throw utils::wecall_exception( "Malformatted bed file line " + std::to_string( lineNumber ) );
        }
    }
}
}