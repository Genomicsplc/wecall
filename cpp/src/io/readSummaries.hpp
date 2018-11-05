// All content Copyright (C) 2018 Genomics plc
#ifndef READ_SUMMARIES_HPP
#define READ_SUMMARIES_HPP

#include "io/read.hpp"
#include "caller/callSet.hpp"
#include "variant/type/variant.hpp"
#include "variant/variantContainer.hpp"
#include "io/readRange.hpp"
#include "utils/interval.hpp"

#include <string>

namespace wecall
{
namespace io
{
    namespace readsummaries
    {
        typedef std::vector< std::vector< int64_t > > coverageDeltas_t;
        struct readCoverage_t
        {
            caller::Region region;
            int64_t minTotalCoverage;
            std::vector< int64_t > minCoveragePerSample;
            std::vector< double_t > averageCoveragePerSample;
        };

        /// Generate a summary of per-base level average read coverage in a given interval
        ///
        /// @param reads Per-sample sets of reads covering the interval
        /// @param beg The start of the interval (included in return value)
        /// @param end The end of the interval (not included in return value, i.e. half-open)

        coverageDeltas_t getReadCoverageDeltas( const perSampleRegionsReads_t & reads,
                                                const caller::Region & subRegion,
                                                const std::size_t nSamples,
                                                const std::size_t nLoci );

        std::vector< readCoverage_t > getChunkedReferenceCalls(
            const caller::Region & subRegion,
            size_t nSamples,
            size_t nLoci,
            const std::vector< std::vector< int64_t > > & readDeltas,
            const double readQualityDeltaThreshold );

        std::vector< readCoverage_t > summariseReadCoverageAndSplitInChunks( const perSampleRegionsReads_t & reads,
                                                                             const caller::Region subRegion,
                                                                             const double readQualityDeltaThreshold );
    }
}
}

#endif
