// All content Copyright (C) 2018 Genomics plc
#ifndef JOB_REDUCE_HPP
#define JOB_REDUCE_HPP

#include <boost/filesystem.hpp>
#include <vector>

#include "utils/timer.hpp"
#include "caller/params.hpp"

namespace echidna
{
namespace caller
{
    class JobReduce
    {
    public:
        JobReduce( const caller::params::Reduce & reduceParams );

        void process();

    private:
        void writeRecords( std::ofstream & out ) const;
        void writeHeader( std::ofstream & out ) const;
        void cleanUp() const;

    private:
        const caller::params::Reduce m_reduceParams;
        std::vector< boost::filesystem::path > m_inputVCFFilePaths;
        utils::timerPtr_t m_timer;
    };
}
}

#endif
