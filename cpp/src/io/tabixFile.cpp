// All content Copyright (C) 2018 Genomics plc
#include <string>
#include "io/tabixFile.hpp"
#include <boost/algorithm/string.hpp>

#include <tabix/tabix.h>

namespace wecall
{
namespace io
{
    TabixFile::TabixFile( std::string filename, std::string indexFilename )
        : m_tabixFile( ti_open( filename.c_str(), indexFilename.c_str() ) )
    {
        this->readHeader();
    }

    TabixFile::~TabixFile() { ti_close( m_tabixFile ); }

    void TabixFile::readHeader()
    {
        const std::size_t N_BYTES = 4096;
        static_assert( N_BYTES > 1, "doesn't make sense to read no data" );
        char data[N_BYTES];
        bgzf_seek( m_tabixFile->fp, 0, SEEK_SET );

        std::string current;
        std::vector< std::string > lines;
        bool complete = false;
        while ( not complete and bgzf_read( m_tabixFile->fp, &data, N_BYTES - 1 ) )
        {
            data[N_BYTES - 1] = '\0';
            current += data;
            if ( current.end() != std::find( current.begin(), current.end(), '\n' ) )
            {
                std::vector< std::string > currentLines;
                boost::algorithm::split( currentLines, current, boost::is_any_of( "\n" ) );
                current = currentLines.back();
                currentLines.pop_back();
                for ( const auto & line : currentLines )
                {
                    if ( boost::starts_with( line, "#" ) )
                    {
                        lines.push_back( line );
                    }
                    else
                    {
                        complete = true;
                        break;
                    }
                }
            }
        }
        m_headerLines = lines;
    }

    std::vector< std::string > TabixFile::fetch( const caller::Region & region ) const
    {
        std::vector< std::string > lines = {};

        if ( ti_iter_t iter = ti_querys( m_tabixFile, region.toString().c_str() ) )
        {
            int len;
            while ( const char * s = ti_read( m_tabixFile, iter, &len ) )
            {
                lines.emplace_back( s );
            }
            ti_iter_destroy( iter );
        }

        return lines;
    }
}
}
