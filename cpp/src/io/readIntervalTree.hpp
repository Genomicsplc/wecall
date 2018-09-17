// All content Copyright (C) 2018 Genomics plc
#ifndef READ_INTERVAL_TREE_HPP
#define READ_INTERVAL_TREE_HPP

#include "io/read.hpp"
#include "utils/intervalTree.hpp"

#include <cstdint>

namespace echidna
{
namespace io
{
    class ReadAlignedEndPosComp;
    class ReadStartPosComp;
    /// interval tree typedefs. Issues when int replaced with int64_t in readIntervalTree_t
    using readIntervalTree_t = utils::IntervalTree< Read, ReadStartPosComp, ReadAlignedEndPosComp >;
    using readIt_t = readIntervalTree_t::iterator;

    /// Utility class for use ordering internal multisets in interval trees.
    /// Note:- Reverse ordering on EndPos
    class ReadAlignedEndPosComp
    {
    public:
        bool operator()( readPtr_t left, readPtr_t right ) const
        {
            return this->operator()( left ) > this->operator()( right );
        }

        int64_t operator()( const readPtr_t & readPtr ) const
        {
            return readPtr->getAlignedEndPos() + readPtr->getLengthAfterAlignedEndPos();
        }
    };

    /// Utility class for use ordering interal multiset for start position in interval tree
    class ReadStartPosComp
    {
    public:
        bool operator()( const readPtr_t & left, const readPtr_t & right ) const
        {
            return this->operator()( left ) < this->operator()( right );
        }

        int64_t operator()( const readPtr_t & readPtr ) const
        {
            return readPtr->getStartPos() - readPtr->getLengthBeforeAlignedStartPos();
        }
    };
}  // io
}  // echidna

#endif  // WECALL_READINTERVALTREE_HPP_H
