// All content Copyright (C) 2018 Genomics plc
#include <chrono>
#include <sstream>

#include "common.hpp"
#include "timer.hpp"
#include "utils/logging.hpp"

namespace echidna
{
namespace utils
{
    std::string encodeString( std::string rawString )
    {
        const std::map< char, std::string > escape{
            {'\\', "\\\\"},
            {'\"', "\\\""},
            {'\a', "\\a"},
            {'\b', "\\b"},
            {'\n', "\\n"},
            {'\t', "\\t"},
            {'\r', "\\r"},
            {'\v', "\\v"},
        };
        std::string encodedString;
        encodedString.push_back( '"' );
        for ( auto c : rawString )
        {
            auto it = escape.find( c );
            if ( it != escape.end() )
            {
                encodedString.append( it->second );
            }
            else
            {
                encodedString.push_back( c );
            }
        }
        encodedString.push_back( '"' );
        return encodedString;
    }

    Timer::Timer( std::string type, std::map< std::string, std::string > metadata )
        : m_type( type ), m_metadata( metadata ), m_duration( 0 )
    {
    }

    void Timer::start() { m_start = std::chrono::steady_clock::now(); }

    void Timer::pause()
    {
        auto end = std::chrono::steady_clock::now();
        m_duration += std::chrono::duration_cast< std::chrono::microseconds >( end - m_start ).count();
    }

    std::string Timer::format_metadata() const
    {
        std::stringstream metadata;
        for ( auto it = m_metadata.cbegin(); it != m_metadata.cend(); )
        {

            metadata << it->first << "=" << encodeString( it->second ) << ";";
            if ( ++it == m_metadata.cend() )
            {
                break;
            }
            metadata << " ";
        }
        return metadata.str();
    }

    Timer::~Timer()
    {
        ECHIDNA_LOG( TIMING, m_type << " " << std::to_string( m_duration ) << "us: " << this->format_metadata() );
    }

    std::map< std::string, std::string > fileMetaData( std::string filename ) { return {{"file", filename}}; }
}
}
