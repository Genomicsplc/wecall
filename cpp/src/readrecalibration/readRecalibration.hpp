// All content Copyright (C) 2018 Genomics plc
#ifndef SLIP_SLIDE_HPP
#define SLIP_SLIDE_HPP

#include "version/version.hpp"
#include "io/bamFile.hpp"
#include "io/fastaFile.hpp"
#include "caller/region.hpp"
#include "io/read.hpp"
#include "io/readRange.hpp"
#include "caller/region.hpp"

#include "readrecalibration/errorCorrectionParameters.hpp"
#include "siteKmerDistribution.hpp"
#include "readrecalibration/kmerDistribution.hpp"
#include "readrecalibration/readDataForErrorPosterior.hpp"
#include "readrecalibration/commonTypes.hpp"

namespace wecall
{
namespace corrector
{
    void floorLowQualityScores( const io::perSampleRegionsReads_t & allReadsInRegion, char qualityFloor, char floorTo );

    void recalibrateDephasingErrors( const io::perSampleRegionsReads_t & allReadsInRegion,
                                     const wecall::io::FastaFile & fa,
                                     const caller::Region & region,
                                     const ErrorCorrectionParameters & errorCorrectionParameters );

    void recalibrateReads( const io::perSampleRegionsReads_t & allReadsInRegion,
                           const wecall::io::FastaFile & fa,
                           const caller::Region & region,
                           const ErrorCorrectionParameters & errorCorrectionParameters );
}
}

#endif
