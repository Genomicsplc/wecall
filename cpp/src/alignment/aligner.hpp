// All content Copyright (C) 2018 Genomics plc
#ifndef ALIGNER_HPP
#define ALIGNER_HPP

#include "common.hpp"

namespace wecall
{
namespace io
{
    class Read;
}  // namespace io

namespace mapping
{
    class HashMapper;
}  // namespace mapping

namespace alignment
{
    class GAlign;

    //-----------------------------------------------------------------------------------------

    double computeLikelihoodForReadAndHaplotype( const io::Read & theRead,
                                                 const int64_t hintPosition,
                                                 const mapping::HashMapper & mapper,
                                                 const alignment::GAlign & aligner );
    //-----------------------------------------------------------------------------------------
}
}

#endif
