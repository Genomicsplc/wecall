// All content Copyright (C) 2018 Genomics plc
#ifndef WECALL_BASE_HPP
#define WECALL_BASE_HPP

#include <boost/program_options.hpp>
#include <string>

namespace wecall
{
class weCallBase
{
public:
    weCallBase( std::string programName )
        : m_programName( programName ),
          m_programVersion( std::string( g_PRODUCT_VERSION ) ),
          m_publicOpts( m_programName + " v" + m_programVersion + " (" + std::string( g_GIT_SHA1, 10 ) +
                        ") configuration parameters" )
    {
    }

    virtual ~weCallBase() {}

protected:
    std::string m_programName;
    std::string m_programVersion;

    boost::program_options::options_description m_publicOpts;
};
}
#endif
