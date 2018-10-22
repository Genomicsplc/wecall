// All content Copyright (C) 2018 Genomics plc
#include "alignment/galign.hpp"

#include <cassert>
#include <algorithm>
#include <iostream>

#include "align.hpp"
#include "utils/sequence.hpp"
#include "utils/logging.hpp"
#include "utils/interval.hpp"

namespace wecall
{
namespace alignment
{
    utils::Interval allowableStartPositionsForAlignment( const int64_t haplotypeLength,
                                                         const int64_t readLength,
                                                         const int paddingLength )
    {
        return utils::Interval( 0, haplotypeLength - readLength + 1 ).getPadded( -paddingLength );
    }

    localGapOpenPenalties_t computeGapOpen( const utils::BasePairSequence & haplotypeSequence,
                                            const errorModel_t & errorModel )
    {
        localGapOpenPenalties_t localGapOpen;
        localGapOpen.resize( haplotypeSequence.size(), 0 );

        // Get gap opening scores along sequence 1, using a homopolymer diploid.
        const auto maxModelIndex = errorModel.size() - 1;

        std::size_t idx = haplotypeSequence.size();
        std::size_t lenHomopolymer = 0;
        char prevBase = constants::gapChar;

        while ( idx-- )
        {
            if ( haplotypeSequence[idx] == prevBase and prevBase != constants::gapChar )
            {
                ++lenHomopolymer;
            }
            else
            {
                lenHomopolymer = 0;
            }

            const auto errorModelIndex = std::min( lenHomopolymer, maxModelIndex );
            const auto gapOpenPenalty = errorModel[errorModelIndex];

            localGapOpen[idx] = gapOpenPenalty;
            prevBase = haplotypeSequence[idx];
        }

        return localGapOpen;
    }

    //-----------------------------------------------------------------------------------------

    GAlign::GAlign( const utils::BasePairSequence & haplotypeSequence,
                    const unsigned short gapExtend,
                    const unsigned short nucleotidePrior,
                    const localGapOpenPenalties_t & localGapOpen )
        : m_haplotypeSequence( haplotypeSequence ),
          m_gapExtend( gapExtend ),
          m_nucPrior( nucleotidePrior ),
          m_localGapOpen( localGapOpen )
    {
        assert( m_localGapOpen.size() == m_haplotypeSequence.size() );
    }

    //-----------------------------------------------------------------------------------------

    int GAlign::computeAlignmentPhredScore( const wecall::utils::BasePairSequence & readSeq,
                                            const wecall::utils::QualitySequence & qual,
                                            const int pos,
                                            char * aln1,
                                            char * aln2 ) const
    {
        assert( pos < 100000 );
        ECHIDNA_ASSERT( qual.size() == readSeq.size(),
                        "Aligner was called with qual string of the wrong length:" + std::to_string( qual.size() ) +
                            " when it should be " + std::to_string( readSeq.size() ) + " to match the read length" );

        const unsigned int readLength = static_cast< int >( readSeq.size() );
        {
            const int haplotypeLength = static_cast< int >( m_haplotypeSequence.size() );
            const auto goodStartPositions =
                allowableStartPositionsForAlignment( haplotypeLength, readLength, constants::needlemanWunschPadding );

            ECHIDNA_ASSERT( goodStartPositions.contains( pos ), "Aligner provided an invalid position to align to: " +
                                                                    std::to_string( pos ) + " not contained in " +
                                                                    goodStartPositions.toString() );
        }

        // the bottom-left and top-right corners of the DP table are just
        // included at the extreme ends of the diagonal, which measures
        // n=8 entries diagonally across.  This fixes the length of the
        // longer (horizontal) sequence to 15 (2*8-1) more than the shorter
        const unsigned int doublePadding = 2 * constants::needlemanWunschPadding - 1;

        const unsigned int haplotypeSegmentLength = readLength + doublePadding;
        const int offset = pos - constants::needlemanWunschPadding;
        auto haplotypeSegmentForAlignment = m_haplotypeSequence.cbegin() + offset;

        int firstPos;

        return needlemanWunschAlignment( haplotypeSegmentForAlignment, readSeq.cbegin(), qual.c_str(),
                                         haplotypeSegmentLength, readLength, m_gapExtend, m_nucPrior, m_localGapOpen,
                                         aln1, aln2, &firstPos, offset );
    }

    //-----------------------------------------------------------------------------------------
}
}
