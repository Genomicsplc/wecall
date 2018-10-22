// All content Copyright (C) 2018 Genomics plc
#include "alignment/aligner.hpp"
#include "io/read.hpp"
#include "mapping/hashMapper.hpp"
#include "alignment/galign.hpp"
#include "common.hpp"

#include <string>
#include <vector>
#include <cmath>
#include <limits>

//_________________________________________________________________________________________________

namespace wecall
{
namespace alignment
{
    double computeLikelihoodForReadAndHaplotype( const io::Read & theRead,
                                                 const int64_t hintPosition,
                                                 const mapping::HashMapper & mapper,
                                                 const alignment::GAlign & aligner )
    {
        /// Compute the best possible alignment score for this read/haplotype pair, and
        /// return it.

        const auto & readSeq = theRead.sequence();
        const auto mapPositions = mapper.mapSequence( readSeq, int64_to_sizet( hintPosition ) );

        if ( mapPositions.empty() )
        {
            return 0.0;
        }
        else
        {
            const auto & readQuals = theRead.getQualities();
            int bestScore = std::numeric_limits< int >::max();

            for ( const auto candidatePos : mapPositions )
            {
                const auto score = aligner.computeAlignmentPhredScore( readSeq, readQuals, candidatePos );
                bestScore = std::min( bestScore, score );
            }

            const auto mapq = theRead.getMappingQuality();
            const double probMappingWrong = stats::fromPhredQ( mapq );
            const double bestAlignmentLikelihood = stats::fromPhredQ( bestScore );
            return bestAlignmentLikelihood * ( 1.0 - probMappingWrong ) + probMappingWrong * 10.0e-20;
        }
    }
}
}

//_________________________________________________________________________________________________
