// All content Copyright (C) 2018 Genomics plc
#ifndef WECALL_MAP_HPP
#define WECALL_MAP_HPP

#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/construct.hpp>

#include "caller/jobReduce.hpp"
#include "common.hpp"
#include "caller/job.hpp"
#include "version/version.hpp"
#include "weCallBase.hpp"

namespace echidna
{
using namespace boost::program_options;
using namespace io;

class weCallMapAndReduce : public weCallBase
{
public:
    weCallMapAndReduce();
    int processJob( int argc, char * argv[] );

private:
    void initOptions();

private:
    options_description m_configOpts;
    options_description m_cmdLineOpts;
};
}

#endif