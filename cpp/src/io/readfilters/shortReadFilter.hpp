// All content Copyright (C) 2018 Genomics plc
#ifndef SHORTREAD_FILTER_HPP
#define SHORTREAD_FILTER_HPP

#include "io/readfilters/readFilter.hpp"

namespace wecall
{
namespace io
{
    /// Simply filters out reads where fragment < 1 read length
    class ShortReadFilter : public ReadFilter
    {
    public:
        /// Constructor
        ///
        ShortReadFilter();
        virtual ~ShortReadFilter() {}
        std::string toString() const override;

    private:
        virtual bool passesFilter_impl( const io::Read & theRead ) override;
    };
}
}

#endif
