// All content Copyright (C) 2018 Genomics plc
#ifndef BOOLEAN_FILTER_HPP
#define BOOLEAN_FILTER_HPP

#include "io/readfilters/readFilter.hpp"
#include <functional>

namespace wecall
{
namespace io
{
    /// Generic filter which takes a boolean function object and evaluates it for each
    /// read.
    class BooleanFilter : public ReadFilter
    {
    public:
        /// Constructor.
        ///
        /// @param theFunction. A function which takes read and returns a boolean
        /// @param filterName The name of this filter
        /// @param negate If true, use the negation of the supplied functions return value to evaluate pass/fail
        BooleanFilter( std::function< const bool(const Read &)> theFunction,
                       const std::string & filterName,
                       const bool negate = false );

        virtual ~BooleanFilter(){};
        std::string toString() const override;

    private:
        virtual bool passesFilter_impl( const io::Read & theRead ) override;

        std::function< const bool(const Read &)> m_function;  /// The function to evaluate
        const std::string m_filterName;  /// The name of this filter (mainly for reporting/debugging)
        const bool m_negate;             /// If true, return the negation of the function
    };
}
}

#endif
