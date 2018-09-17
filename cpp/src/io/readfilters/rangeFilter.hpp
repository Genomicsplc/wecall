// All content Copyright (C) 2018 Genomics plc
#ifndef RANGE_FILTER_HPP
#define RANGE_FILTER_HPP

#include "io/readfilters/readFilter.hpp"
#include <functional>
#include <sstream>

namespace echidna
{
namespace io
{
    /// Generic filter which takes a function object and evaluates it for each
    /// read, returning true if the result falls within the specified range.
    template < typename T >
    class RangeFilter : public ReadFilter
    {
    public:
        /// Constructor
        ///
        /// @param theFunction A callable object which takes a read and returns a T
        /// @param filterName The name of this filter
        /// @param min The minimum acceptable value
        /// @param max The maximum acceptable value
        RangeFilter( std::function< const T(const Read &)> theFunction,
                     const std::string & filterName,
                     const T min,
                     const T max )
            : m_function( theFunction ), m_filterName( filterName ), m_min( min ), m_max( max )
        {
        }

        virtual ~RangeFilter() {}

        std::string toString() const
        {
            std::stringstream repr;
            repr << this->m_filterName << "(min = " << m_min << ", max = " << m_max << ")";
            return repr.str();
        }

    private:
        virtual bool passesFilter_impl( const io::Read & theRead )
        {
            const T result = this->m_function( theRead );

            if ( result >= this->m_min && result <= this->m_max )
            {
                return true;
            }

            else
            {
                return false;
            }
        }

        std::function< const T(const Read &)> m_function;  /// The function to evaluate
        const std::string m_filterName;                    /// The name of this filter (mainly for reporting/debugging)
        const T m_min;                                     /// Lower bound
        const T m_max;                                     /// Upper bound
    };
}
}

#endif
