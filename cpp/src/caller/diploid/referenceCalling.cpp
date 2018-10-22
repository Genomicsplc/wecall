// All content Copyright (C) 2018 Genomics plc
#include "referenceCalling.hpp"

#include "common.hpp"
#include "io/readSummaries.hpp"
#include "stats/models.hpp"

namespace wecall
{
namespace caller
{
    namespace model
    {
        //-----------------------------------------------------------------------------------------

        std::vector< Call > buildRefCall( caller::Region refInterval,
                                          const io::perSampleRegionsReads_t & reads,
                                          double maxUncalledVarQ,
                                          const std::vector< std::size_t > & ploidyPerSample,
                                          const double readQualityDeltaThreshold )
        {
            // First get coverage data for each sample

            ECHIDNA_ASSERT( refInterval.size() > 0, "Can not call reference on empty reference interval." );

            const auto nSamples = reads.size();
            const auto readCoverageChunks = io::readsummaries::summariseReadCoverageAndSplitInChunks(
                reads, refInterval, readQualityDeltaThreshold );

            // Calculate quality of reference call as follows:
            //
            // 1.    If no post read-filtering evidence of variants in the region, calculate the probability that
            //        we may only have sequenced one copy of a diploid, thus possibly missing a variant
            //        then convert to a PHRED quality score.
            // 2.    If some evidence of variants, but variants rejected during cluster processing, use the
            //        maximum probability of the rejected vairants to calculate the PHRED quality score and then
            //        take the minimum of the two quality scores.

            std::vector< Call > calls;
            for ( auto const & chunk : readCoverageChunks )
            {
                double refQ = getRefQFromCoverageRow( chunk.minTotalCoverage );
                if ( not std::isnan( refQ ) and maxUncalledVarQ > 0.0 )
                {
                    double notVarQ = stats::toPhredQ( 1.0 - stats::fromPhredQ( maxUncalledVarQ ) );
                    refQ = std::min( refQ, notVarQ );
                }

                genoCalls_t genoCalls;
                for ( const auto & ploidy : ploidyPerSample )
                {
                    genoCalls.emplace_back( ploidy, caller::Call::REF );
                }

                // Create a reference call
                Call call( nullptr, chunk.region.interval(), refQ, nSamples, genoCalls );

                // Add annotations - note: start positions is inclusive, end is exclusive and both
                // are 0-indexed, but VCF requires an inclusive 1-indexed range.
                call.addAnnotation( Annotation::BEG, chunk.region.start() + 1 );  // make it 1-based
                call.addAnnotation( Annotation::END, chunk.region.end() );  // make it 1-based with inclusive interval
                call.addAnnotation( Annotation::LEN, chunk.region.size() );

                // And the sample specific annotations (note: sample coverage starts at row 1)
                for ( std::size_t i = 0; i < nSamples; ++i )
                {
                    call.samples[i].addAnnotation(
                        Annotation::FORMAT_DP, static_cast< int64_t >( round( chunk.averageCoveragePerSample[i] ) ) );
                    call.samples[i].addAnnotation( Annotation::MIN_DP, chunk.minCoveragePerSample[i] );
                    call.samples[i].addAnnotation( Annotation::GQ,
                                                   getRefQFromCoverageRow( chunk.minCoveragePerSample[i] ) );
                }
                calls.push_back( call );
            }
            return calls;
        }

        //-----------------------------------------------------------------------------------------

        double getQualityFromCoverage( const std::int64_t minCoverage )
        {
            return stats::toPhredQ( stats::probAllReadsFromSingleCopyOfDiploid( static_cast< int >( minCoverage ) ) );
        }

        //-----------------------------------------------------------------------------------------

        double getRefQFromCoverageRow( const std::int64_t minCoverage )
        {
            if ( minCoverage == 0 )
            {
                // Special case (to avoid returning something which is fractionally negative)
                return std::numeric_limits< double >::quiet_NaN();
            }
            else
            {
                return getQualityFromCoverage( minCoverage );
            }
        }
    }
}
}
