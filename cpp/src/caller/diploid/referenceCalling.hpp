// All content Copyright (C) 2018 Genomics plc
#ifndef WECALL_REFERENCECALLING_H
#define WECALL_REFERENCECALLING_H

#include "caller/callSet.hpp"
#include "io/readRange.hpp"

namespace echidna
{
namespace caller
{
    namespace model
    {
        std::vector< Call > buildRefCall( caller::Region refInterval,
                                          const io::perSampleRegionsReads_t & reads,
                                          double maxUncalledVarQ,
                                          const std::vector< std::size_t > & ploidy,
                                          const double readQualityDeltaThreshold );

        double getQualityFromCoverage( const int64_t minCoverage );
        double getRefQFromCoverageRow( const int64_t minCoverage );
    }
}
}

#endif  // WECALL_REFERENCECALLING_H
