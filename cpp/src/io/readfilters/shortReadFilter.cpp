// All content Copyright (C) 2018 Genomics plc
#include "io/readfilters/shortReadFilter.hpp"
#include <sstream>
#include <cstdlib>

namespace echidna
{
namespace io
{
    //-----------------------------------------------------------------------------------------

    ShortReadFilter::ShortReadFilter() {}

    //-----------------------------------------------------------------------------------------

    bool ShortReadFilter::passesFilter_impl( const io::Read & theRead )
    {
        const auto insertSize = std::abs( theRead.getInsertSize() );
        const auto readLengthNoSoftClipping = theRead.cigar().lengthInSeqWithoutSoftClipping();

        return insertSize >= readLengthNoSoftClipping;
    }

    //-----------------------------------------------------------------------------------------

    std::string ShortReadFilter::toString() const
    {
        std::stringstream repr;
        repr << "ShortReadFilter";
        return repr.str();
    }

    //-----------------------------------------------------------------------------------------
}
}

// std::function<const int(const Read&)> f1 = &Read::getInsertSize;
// std::function<const int(const Read&)> f2 = &Read::getLength;
// std::function<const bool(const Read& theRead)> f3 = [f1, f2](const Read& theRead){ return f1(theRead) < f2(theRead);
// };
