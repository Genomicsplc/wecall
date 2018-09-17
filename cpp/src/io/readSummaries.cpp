// All content Copyright (C) 2018 Genomics plc
#include "io/readSummaries.hpp"
#include "caller/diploid/referenceCalling.hpp"
#include <numeric>

namespace echidna
{
namespace io
{
    namespace readsummaries
    {

        //-------------------------------------------------------------------------------------

        coverageDeltas_t getReadCoverageDeltas( const perSampleRegionsReads_t & reads,
                                                const caller::Region & subRegion,
                                                const std::size_t nSamples,
                                                const std::size_t nLoci )
        {
            // Record starting read depath and changes in read depth in region.
            coverageDeltas_t readDeltas( nSamples, std::vector< int64_t >( nLoci, 0 ) );

            std::size_t sampleIndex = 0;
            for ( const auto & sampleReads : reads )
            {
                const auto restrictedRegion = sampleReads.second.getRegions().getSpan().getIntersect( subRegion );
                const auto regionsReads = sampleReads.second.getSubRegionReads( restrictedRegion );

                for ( const auto & read : regionsReads )
                {
                    auto readBeg = read.getStartPos();
                    auto readEnd = read.getAlignedEndPos();

                    // Because reads cannot be sorted by end as well as start, we may have some
                    // near the beginning of the range which fall short after alignment - ignore them.
                    if ( readEnd <= subRegion.start() )
                    {
                        continue;
                    }

                    //  Accumulate coverage at first locus
                    if ( readBeg <= subRegion.start() )
                    {
                        ++readDeltas[sampleIndex][0];
                    }

                    // Record deltas where a read stops or starts
                    if ( readBeg > subRegion.start() and readBeg < subRegion.end() )
                    {
                        ++readDeltas[sampleIndex][readBeg - subRegion.start()];
                    }
                    if ( readEnd < subRegion.end() )
                    {
                        --readDeltas[sampleIndex][readEnd - subRegion.start()];
                    }
                }
                sampleIndex++;
            }
            return readDeltas;
        }

        //-------------------------------------------------------------------------------------

        std::vector< readCoverage_t > getChunkedReferenceCalls( const caller::Region & subRegion,
                                                                const size_t nSamples,
                                                                const size_t nLoci,
                                                                const coverageDeltas_t & readDeltas,
                                                                const double readQualityDeltaThreshold )
        {
            std::vector< readCoverage_t > readCoverageChunks;
            std::vector< int64_t > totalCoverageOverSamples( std::vector< int64_t >( nLoci, 0 ) );

            std::vector< int64_t > currentReadDepthPerSample;
            for ( size_t sampleIndex = 0; sampleIndex < nSamples; sampleIndex++ )
            {
                currentReadDepthPerSample.push_back( readDeltas[sampleIndex][0] );
            }
            std::vector< int64_t > minCoveragePerSample = currentReadDepthPerSample;
            std::vector< int64_t > totalReadsAccPerSample = currentReadDepthPerSample;
            std::vector< double_t > averageReadDepthInChunkPerSample( std::vector< double_t >( nSamples, 0 ) );

            totalCoverageOverSamples[0] =
                std::accumulate( currentReadDepthPerSample.begin(), currentReadDepthPerSample.end(), 0 );

            int64_t numPositions = 1;
            int64_t startPosChunk = 0;

            // create chunks of reference calls
            for ( std::size_t i = 1; i < nLoci; ++i )
            {

                // decide if we need to start a new chunk
                bool startNewChunk = false;
                for ( size_t sampleIndex = 0; sampleIndex < nSamples; sampleIndex++ )
                {
                    averageReadDepthInChunkPerSample[sampleIndex] =
                        totalReadsAccPerSample[sampleIndex] / static_cast< float_t >( numPositions );

                    // advance the current read depth by one position
                    currentReadDepthPerSample[sampleIndex] =
                        currentReadDepthPerSample[sampleIndex] + readDeltas[sampleIndex][i];

                    // start new chunk if read depth changes significantly
                    double minQuality = caller::model::getQualityFromCoverage( minCoveragePerSample[sampleIndex] );
                    double currentQuality =
                        caller::model::getQualityFromCoverage( currentReadDepthPerSample[sampleIndex] );
                    double relativeQualityChange = ( currentQuality - minQuality ) / minQuality;

                    // at the moment the quality calculation is instable for high read counts JPLAT-933
                    if ( std::abs( relativeQualityChange ) > readQualityDeltaThreshold )
                    {
                        startNewChunk = true;
                    }
                }

                // create a new chunk if there is a strong change in read depth or for last entry
                if ( startNewChunk )
                {
                    // finish creation of read chunk
                    int64_t minTotalCoverage =
                        *std::min_element( minCoveragePerSample.begin(), minCoveragePerSample.end() );
                    readCoverage_t readCoverage = {
                        caller::Region( subRegion.contig(), subRegion.start() + startPosChunk,
                                        subRegion.start() + i ),  // half-open interval
                        minTotalCoverage,
                        minCoveragePerSample,
                        averageReadDepthInChunkPerSample};
                    readCoverageChunks.push_back( readCoverage );

                    // reset accumulators
                    minCoveragePerSample = currentReadDepthPerSample;
                    totalReadsAccPerSample = currentReadDepthPerSample;
                    numPositions = 1;
                    startPosChunk = i;
                }
                else
                {
                    // move one position along region
                    for ( size_t sampleIndex = 0; sampleIndex < nSamples; sampleIndex++ )
                    {
                        totalCoverageOverSamples[i] += currentReadDepthPerSample[sampleIndex];
                        totalReadsAccPerSample[sampleIndex] += currentReadDepthPerSample[sampleIndex];
                        minCoveragePerSample[sampleIndex] =
                            std::min( currentReadDepthPerSample[sampleIndex], minCoveragePerSample[sampleIndex] );
                    }
                    ++numPositions;
                }
            }

            // add the last chunk
            int64_t minTotalCoverage = *std::min_element( minCoveragePerSample.begin(), minCoveragePerSample.end() );
            for ( size_t sampleIndex = 0; sampleIndex < nSamples; sampleIndex++ )
            {
                averageReadDepthInChunkPerSample[sampleIndex] =
                    totalReadsAccPerSample[sampleIndex] / static_cast< float_t >( numPositions );
            }
            readCoverage_t readCoverage = {
                caller::Region( subRegion.contig(), subRegion.start() + startPosChunk, subRegion.end() ),
                minTotalCoverage,
                minCoveragePerSample,
                averageReadDepthInChunkPerSample};
            readCoverageChunks.push_back( readCoverage );

            return readCoverageChunks;
        }

        //-------------------------------------------------------------------------------------

        std::vector< readCoverage_t > summariseReadCoverageAndSplitInChunks( const perSampleRegionsReads_t & reads,
                                                                             const caller::Region subRegion,
                                                                             const double qualityDeltaThreshold )
        {
            std::size_t nSamples = reads.size();
            std::size_t nLoci = int64_to_sizet( subRegion.size() );

            // Initialise coverage structure
            // read deltas is holding in the first position the read depth at this position and then the delta of reads
            // in each subsequent field
            std::vector< std::vector< int64_t > > readDeltas =
                getReadCoverageDeltas( reads, subRegion, nSamples, nLoci );

            std::vector< readCoverage_t > readCoverageChunks =
                getChunkedReferenceCalls( subRegion, nSamples, nLoci, readDeltas, qualityDeltaThreshold );

            return readCoverageChunks;
        }
    }
}
}
