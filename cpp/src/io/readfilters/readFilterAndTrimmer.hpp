// All content Copyright (C) 2018 Genomics plc
#ifndef READ_FILTER_AND_TRIMMER_HPP
#define READ_FILTER_AND_TRIMMER_HPP

#include <vector>

#include "io/read.hpp"
#include "common.hpp"
#include "caller/params.hpp"
#include "io/readfilters/readFilter.hpp"

namespace wecall
{
namespace io
{
    /// Manages all read filters and trimmers that are applied to individual reads during reading
    class ReadFilterAndTrimmer
    {
    public:
        ReadFilterAndTrimmer( const caller::params::Filters & filterParams );

        /// Trim the read and return true if passed all filters and the trimmed read has length > 0
        bool trimAndFilter( readPtr_t read ) const;

    private:
        void trim( readPtr_t read ) const;
        bool passesFilters( readPtr_t read ) const;
        bool hasLength( readPtr_t read ) const;
        bool isSimilarToPrevious( readPtr_t read ) const;

        std::vector< ReadFilterPtr_t > m_filters;

        bool m_overlapTrim;
        bool m_shortReadTrim;
        bool m_noSimilarReads;

        mutable readPtr_t m_previousRead = nullptr;
    };
}
}

#endif  // READ_FILTERS_MANAGER_HPP
