// All content Copyright (C) 2018 Genomics plc
#include <iterator>
#include <algorithm>
#include <fstream>

#include <boost/filesystem.hpp>

#include "utils/timer.hpp"
#include "vcf/header.hpp"
#include "caller/jobReduce.hpp"
#include "utils/logging.hpp"
#include "caller/params.hpp"

namespace wecall
{
namespace caller
{
    namespace fs = boost::filesystem;

    JobReduce::JobReduce( const caller::params::Reduce & reduceParams )
        : m_reduceParams( reduceParams ),
          m_timer( std::make_shared< utils::Timer >( "IO", utils::fileMetaData( reduceParams.outputDataSink() ) ) )
    {
        WECALL_ERROR( fs::is_directory( m_reduceParams.inputDir() ),
                       m_reduceParams.inputDir() + " is not a directory" );

        for ( fs::directory_iterator directory_iterator( m_reduceParams.inputDir() );
              directory_iterator != fs::directory_iterator(); ++directory_iterator )
        {
            fs::path filePath = directory_iterator->path();
            WECALL_ERROR( ( fs::is_regular_file( filePath ) ), filePath.string() + " is not a file" );
            WECALL_ERROR( ( filePath.filename().extension() == fs::path( ".vcf" ) ),
                           "file " + filePath.string() + " is not a VCF" );

            m_inputVCFFilePaths.emplace_back( filePath );
        }

        WECALL_ERROR( ( not m_inputVCFFilePaths.empty() ), "directory " + m_reduceParams.inputDir() + " is empty" );

        std::sort( m_inputVCFFilePaths.begin(), m_inputVCFFilePaths.end() );
    }

    void JobReduce::process()
    {
        utils::ScopedTimerTrigger scopedTimerTrigger( m_timer );
        std::ofstream out( m_reduceParams.outputDataSink(), std::ios_base::out );
        writeHeader( out );
        writeRecords( out );
        out.close();

        cleanUp();

        WECALL_LOG( INFO, "Reduce finished" );
    }

    void JobReduce::cleanUp() const { boost::filesystem::remove_all( m_reduceParams.inputDir() ); }

    void JobReduce::writeHeader( std::ofstream & out ) const
    {
        bool commonContentWritten = false;
        std::string chromLine;

        for ( const auto & file : m_inputVCFFilePaths )
        {
            std::ifstream in( file.string(), std::ios_base::in );

            std::string temp;
            while ( std::getline( in, temp ) )
            {
                if ( temp.substr( 0, vcf::Header::chromKey().size() ) == vcf::Header::chromKey() )
                {
                    chromLine = temp;
                    break;
                }

                if ( not commonContentWritten )
                {
                    out << temp << "\n";
                }
                else if ( temp.substr( 0, vcf::Header::contigKey().size() ) == vcf::Header::contigKey() )
                {
                    out << temp << "\n";
                }

                temp.clear();
            }

            commonContentWritten = true;
        }

        out << chromLine << "\n";
    }

    void JobReduce::writeRecords( std::ofstream & out ) const
    {
        for ( const auto & file : m_inputVCFFilePaths )
        {
            WECALL_LOG( INFO, "Processing " << file );

            std::ifstream in( file.string(), std::ios_base::in );

            std::string temp;

            // look for the line that starts with "#CHROM"
            bool isNextLineARecord = false;
            while ( not isNextLineARecord and std::getline( in, temp ) )
            {
                if ( temp.substr( 0, vcf::Header::chromKey().size() ) == vcf::Header::chromKey() )
                {
                    isNextLineARecord = true;
                }

                temp.clear();
            }

            WECALL_ERROR( isNextLineARecord, "file " + file.string() + " is not a valid VCF" );

            std::istreambuf_iterator< char > beginSource( in );
            std::istreambuf_iterator< char > endSource;
            std::ostreambuf_iterator< char > beginDest( out );
            std::copy( beginSource, endSource, beginDest );

            WECALL_ASSERT( in.good(), "Input data stream in error state" );
            in.peek();
            WECALL_ASSERT( in.eof(), "Input data stream not depleted" );

            in.close();
        }
    }
}
}
