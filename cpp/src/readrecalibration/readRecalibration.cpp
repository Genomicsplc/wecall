// All content Copyright (C) 2018 Genomics plc
#include <string>
#include <cstdint>

#include "readrecalibration/readRecalibration.hpp"
#include "readrecalibration/errorCorrectionParameters.hpp"
#include "io/readRange.hpp"
#include "caller/region.hpp"

namespace echidna
{
namespace corrector
{

    void floorLowQualityScores( const io::perSampleRegionsReads_t & allReadsInRegion, char qualityFloor, char floorTo )
    {
        for ( auto & nameAndReads : allReadsInRegion )
        {
            for ( auto & read : nameAndReads.second )
            {
                for ( auto & qualityCharacter : read.qualities() )
                {
                    if ( qualityCharacter <= qualityFloor )
                    {
                        qualityCharacter = floorTo;
                    }
                }
            }
        }
    }

    void recalibrateDephasingErrors( const io::perSampleRegionsReads_t & allReadsInRegion,
                                     const echidna::io::FastaFile & fa,
                                     const caller::Region & region,
                                     const ErrorCorrectionParameters & errorCorrectionParameters )
    {
        for ( auto & nameAndReads : allReadsInRegion )
        {
            const auto nReads = std::distance( nameAndReads.second.begin(), nameAndReads.second.end() );

            if ( nReads == 0 )
            {
                continue;
            }

            const auto readsStartEndPos = io::readsAlignedStartEnd( nameAndReads.second );

            KmerDistribution kmerDistribution( region.contig(), fa, readsStartEndPos.first, readsStartEndPos.second );

            std::vector< ReadDataForErrorPosterior > dataForErrorPosterior;

            for ( auto read_it = nameAndReads.second.begin(); read_it != nameAndReads.second.end(); ++read_it )
            {
                // Ugly way of retrieving shared_ptr<T> from interval tree. Means only one loop of tree is needed.
                io::readPtr_t readPtr = read_it.getSharedPtr();
                kmerDistribution.updateKmerHistogram( readPtr );
                dataForErrorPosterior.emplace_back( readPtr );
            }

            kmerDistribution.finalise( errorCorrectionParameters.pref );

            for ( auto & readDataForErrorPosterior : dataForErrorPosterior )
            {
                readDataForErrorPosterior.calculateProbabilities( kmerDistribution, errorCorrectionParameters );
            }

            // run HMM across all reads for the first time and accumulate
            kmerDistribution.resetErrorCountData( errorCorrectionParameters.perr );  // populate with pseudo-counts
            kmerDistribution.updateErrorPosteriors();  // this sets error transition probs to an initial value

            for ( auto & readDataForErrorPosterior : dataForErrorPosterior )
            {
                readDataForErrorPosterior.runHmm( kmerDistribution );
            }
            // Use the updated error counts from HMM to update the error posteriors.
            kmerDistribution.updateErrorPosteriors();

            // run HMM across reads for 2nd time; update the per-site posterior probabilities of transitioning into the
            // error state
            // Reset error correction parameters to re-run HMM.
            kmerDistribution.resetErrorCountData( errorCorrectionParameters.perr );
            for ( auto & readDataForErrorPosterior : dataForErrorPosterior )
            {
                readDataForErrorPosterior.runHmm( kmerDistribution );
            }

            for ( auto & readDataForErrorPosterior : dataForErrorPosterior )
            {
                readDataForErrorPosterior.recalibrateRead();
            }
        }
    }

    void recalibrateReads( const io::perSampleRegionsReads_t & allReadsInRegion,
                           const echidna::io::FastaFile & fa,
                           const caller::Region & region,
                           const ErrorCorrectionParameters & errorCorrectionParameters )
    {
        recalibrateDephasingErrors( allReadsInRegion, fa, region, errorCorrectionParameters );
        floorLowQualityScores( allReadsInRegion, errorCorrectionParameters.qualityFloor,
                               errorCorrectionParameters.floorTo );
    }
}
}
