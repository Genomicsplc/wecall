// All content Copyright (C) 2018 Genomics plc
#ifndef GALIGN_HPP
#define GALIGN_HPP

#include <string>
#include <vector>
#include "utils/interval.hpp"
#include "utils/sequence.hpp"
#include "align.hpp"

namespace echidna
{
namespace alignment
{
    localGapOpenPenalties_t computeGapOpen( const utils::BasePairSequence & haplotypeSequence,
                                            const errorModel_t & errorModel );

    utils::Interval allowableStartPositionsForAlignment( const int64_t haplotypeLength,
                                                         const int64_t readLength,
                                                         const int paddingLength );

    class GAlign
    {
    public:
        GAlign( const utils::BasePairSequence & haplotypeSequence,
                const short gapExtend,
                const short nucleotidePrior,
                const localGapOpenPenalties_t & localGapOpen );

        int computeAlignmentPhredScore( const utils::BasePairSequence & readSeq,
                                        const utils::QualitySequence & qual,
                                        const int pos,
                                        char * aln1 = NULL,
                                        char * aln2 = NULL ) const;

    private:
        const utils::BasePairSequence m_haplotypeSequence;
        const short m_gapExtend;
        const short m_nucPrior;
        const localGapOpenPenalties_t m_localGapOpen;
    };
}
}

#endif
