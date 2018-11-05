// All content Copyright (C) 2018 Genomics plc
#ifndef READ_FILTER_HPP
#define READ_FILTER_HPP

#include "common.hpp"
#include "io/read.hpp"

#include <memory>

namespace wecall
{
namespace io
{
    /// Base class for read filters.
    class ReadFilter
    {
    public:
        /// Constructor
        ReadFilter() {}

        /// Destructor
        virtual ~ReadFilter() {}

        /// Returns a string representation of the filter instance
        virtual std::string toString() const = 0;

        /// Return true if the read passes this filter. Otherwise false.
        ///
        /// @param A single read to check
        /// @return The result
        bool passesFilter( const io::Read & theRead );

    private:
        /// Virtual function that implements the filter
        virtual bool passesFilter_impl( const io::Read & theRead ) = 0;
    };

    using ReadFilterPtr_t = std::shared_ptr< ReadFilter >;
}
}

#endif
