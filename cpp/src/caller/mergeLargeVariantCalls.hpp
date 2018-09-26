#ifndef WECALL_MERGELARGEVARIANTCALLS_H
#define WECALL_MERGELARGEVARIANTCALLS_H

#include "caller/callSet.hpp"

namespace echidna
{
namespace caller
{

    class lvcMerger
    {
    private:
    public:
        callVector_t removeReferenceCalls( const callVector_t & calls,
                                           const std::vector< std::size_t > & thisAreaPloidies,
                                           const std::vector< std::size_t > & defaultPloidies ) const;

        std::vector< size_t > getOverlappingLvcCallIndices( const Call & call, const callVector_t & lvcCalls ) const;

        void mergeAndCorrectGenotypes( callVector_t & otherCalls,
                                       callVector_t & lvcCalls,
                                       const std::vector< std::size_t > & thisAreaPloidies,
                                       const std::vector< std::size_t > & defaultPloidy ) const;
    };
}
}

#endif