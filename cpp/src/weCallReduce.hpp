// All content Copyright (C) 2018 Genomics plc
#ifndef WECALL_REDUCE_HPP
#define WECALL_REDUCE_HPP

#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>

#include "caller/jobReduce.hpp"
#include "common.hpp"
#include "caller/job.hpp"
#include "version/version.hpp"
#include "weCallBase.hpp"

namespace echidna
{
class weCallReduce : public weCallBase
{
public:
    weCallReduce();

    int processJob( int argc, char * argv[] );

private:
    void initOptions();
};
}

#endif