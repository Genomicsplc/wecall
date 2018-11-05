// All content Copyright (C) 2018 Genomics plc
#include "io/bamFile.hpp"
#include "utils/exceptions.hpp"

#include <ctype.h>
#include <assert.h>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <memory>
#include <set>
#include <utility>

#include "caller/region.hpp"

#include "boost/algorithm/string.hpp"
#include "boost/filesystem.hpp"
#include "utils/timer.hpp"

namespace wecall
{
namespace io
{

    //-----------------------------------------------------------------------------------------

    BamFile::BamFile( std::string fileName )
        : m_fileName( fileName ),
          m_samFile( nullptr ),
          m_index( nullptr ),
          m_timer( std::make_shared< utils::Timer >( "IO", utils::fileMetaData( fileName ) ) )

    {
        // Dummy commit
        m_index = bam_index_load( fileName.c_str() );

        if ( m_index == nullptr )
        {
            WECALL_LOG( FATAL, "Could not load index for BAM file "
                                    << fileName << ". Check that index file exists and is readable" );
            throw utils::wecall_exception( "Cannot load BAM index" );
        }

        m_samFile = samopen( fileName.c_str(), "rb", nullptr );
        m_samplesByID = this->getSamplesByID();
    }

    //-----------------------------------------------------------------------------------------

    BamFile::~BamFile()
    {
        if ( m_samFile != nullptr )
        {
            samclose( m_samFile );
            m_samFile = nullptr;
        }

        if ( m_index != nullptr )
        {
            bam_index_destroy( m_index );
            m_index = nullptr;
        }
    }

    //-----------------------------------------------------------------------------------------

    std::string BamFile::getContigName( const std::size_t tid ) const
    {
        if ( tid < std::size_t( m_samFile->header->n_targets ) )
        {
            return std::string( m_samFile->header->target_name[tid] );
        }
        else
        {
            throw utils::wecall_exception( "Invalid tid value in BamFile::getContigName" );
        }
    }

    //-----------------------------------------------------------------------------------------

    std::vector< std::string > BamFile::getSampleNames() const
    {
        if ( m_samplesByID.empty() )
        {
            return {this->default_sample_name()};
        }
        else
        {
            std::set< std::string > sampleNames;

            for ( auto & keyVal : m_samplesByID )
            {
                sampleNames.insert( keyVal.second );
            }

            std::vector< std::string > mergedSampleNames( sampleNames.size() );
            copy( sampleNames.begin(), sampleNames.end(), mergedSampleNames.begin() );

            return mergedSampleNames;
        }
    }

    //-----------------------------------------------------------------------------------------

    std::map< std::string, std::string > BamFile::getSamplesByID() const
    {
        const std::string theHeader( m_samFile->header->text );

        std::vector< std::string > headerLines;
        boost::algorithm::split( headerLines, theHeader, boost::is_any_of( "\n" ) );
        std::map< std::string, std::string > samplesByID;

        for ( const auto line : headerLines )
        {
            if ( boost::algorithm::starts_with( line, "@RG" ) )
            {
                std::vector< std::string > columns;
                boost::algorithm::split( columns, line, boost::is_any_of( "\t" ) );

                std::string theID;
                std::string sampleName;

                for ( const auto column : columns )
                {
                    if ( boost::algorithm::starts_with( column, "ID:" ) )
                    {
                        theID = column.substr( 3 );
                    }

                    if ( boost::algorithm::starts_with( column, "SM:" ) )
                    {
                        sampleName = column.substr( 3 );
                    }
                }

                if ( theID.empty() )
                {
                    continue;
                }

                else if ( sampleName.empty() )
                {
                    WECALL_LOG( WARNING, "Could not find sample name for read group with id " << theID );
                    samplesByID[theID] = "UNKNOWN_SAMPLE";
                    continue;
                }

                else
                {
                    WECALL_LOG( DEBUG, "Found sample name " << sampleName << " for read group id " << theID );
                    samplesByID[theID] = sampleName;
                }
            }
        }

        return samplesByID;
    }

    std::string BamFile::default_sample_name() const
    {
        boost::filesystem::path path = boost::filesystem::path( m_fileName );
        boost::filesystem::path fileName = path.filename();
        boost::filesystem::path ext = path.extension();
        std::string sampleName = std::string( fileName.string(), 0, fileName.string().size() - ext.string().size() );
        return sampleName;
    }

    //-----------------------------------------------------------------------------------------

    bamFileIteratorPtr_t BamFile::readRegion( const caller::Region & blockRegion,
                                              utils::referenceSequencePtr_t m_refSequence )
    {
        if ( not m_samFile )
        {
            throw utils::wecall_exception( "Attempted to call 'fetch' on a closed BAM file" );
        }

        int rtid = 0;
        int rstart = 0;
        int rend = 0;

        const caller::Region clippedRegion( blockRegion.contig(), std::max( blockRegion.start(), int64_t( 0L ) ),
                                            blockRegion.end() );
        bam_parse_region( m_samFile->header, clippedRegion.toString().c_str(), &rtid, &rstart, &rend );

        if ( rtid < 0 )
        {
            WECALL_LOG( WARNING, "Attempted to load an invalid contig \"" + clippedRegion.contig() +
                                      "\" from the BAM file - " +
                                      "Check that the contig names in the reference file match those in the BAM." );
            return nullptr;
        }

        if ( rstart > rend )
        {
            WECALL_LOG( FATAL, "Could not retrieve region " << clippedRegion << " from BAM file " << m_fileName );
            WECALL_LOG( FATAL, "Check that the contig range in the FASTA file matches that in the BAM" );
            throw utils::wecall_exception( "Invalid region - BAM file start coordinate > end" );
        }

        bam_fetch_iterator_t * bam_fetch_iterator =
            bam_init_fetch_iterator( m_samFile->x.bam, m_index, rtid, rstart, rend );
        if ( m_samplesByID.empty() )
        {
            return std::make_shared< BamFileWithoutReadGroupIterator >( bam_fetch_iterator, this->default_sample_name(),
                                                                        m_refSequence, m_timer );
        }
        else
        {
            return std::make_shared< BamFileIterator >( bam_fetch_iterator, m_samplesByID, m_refSequence, m_timer );
        }
    }
}
}
