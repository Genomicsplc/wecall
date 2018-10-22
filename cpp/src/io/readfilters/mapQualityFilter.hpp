// All content Copyright (C) 2018 Genomics plc
#ifndef MAPQUAL_FILTER_HPP
#define MAPQUAL_FILTER_HPP

#include "io/readfilters/readFilter.hpp"

namespace wecall
{
namespace io
{
    /// Filter which simply checks the mapping quality of the specified read.
    class MapQualFilter : public ReadFilter
    {
    public:
        /// Constructor.
        ///
        /// @param threshold the minimum acceptable value for mapping quality.
        explicit MapQualFilter( int64_t threshold ) : m_threshold( threshold ) {}

        /// Destructor
        virtual ~MapQualFilter() {}

        /// Convert to string representation
        virtual std::string toString() const override;

    private:
        virtual bool passesFilter_impl( const io::Read & theRead ) override;
        int64_t m_threshold;  /// Min value for mapping quality.
    };
}
}

#endif
