// All content Copyright (C) 2018 Genomics plc
#ifndef BASIC_HASH_HPP
#define BASIC_HASH_HPP

#include <string>
#include <vector>
#include <set>
#include <utils/interval.hpp>
#include "utils/sequence.hpp"

#include "utils/logging.hpp"
#include "boost/optional.hpp"

namespace echidna
{
namespace mapping
{
    namespace constants
    {
        constexpr static int endOfChain = -1;
        constexpr static int repeatChain = -2;
        constexpr static int maxRepeatCount = 10;
        constexpr static double fractionOfKmersToConsider = 1.0 - 1.0 / maxRepeatCount;
    }

    class HashFunction
    {
    public:
        HashFunction( const utils::BasePairSequence & sequence, unsigned int kmerSize );
        int next( const char newBase );

    private:
        const int m_bitShift;
        static constexpr int m_bitMask = 6;  // a,A->000 (0)  c,C->010 (2)  t,T->100 (4)  g,G->110 (6)
        static constexpr int m_nBits = 2;

        int m_currentHash;
    };

    class KmerMatches
    {
    public:
        KmerMatches( const utils::BasePairSequence & paddedHaplotypeSequence,
                     const int kmerSize,
                     const int haplotypePadding );

        std::vector< std::size_t > countKmerMatches( const echidna::utils::BasePairSequence & readSequence ) const;

        utils::Interval allowableStartPositionsInHaplotype( const int64_t readSeqLen ) const;

    private:
        void indexSequence( const utils::BasePairSequence & haplotypeSequence );

        bool isHashHit( const int hashTableValue ) const { return hashTableValue >= 0; }

    private:
        const unsigned int m_paddedHaplotypeSequenceLength;

        const unsigned int m_kmerSize;
        const int m_haplotypePadding;
        const int m_hashSize;

        std::vector< int > m_hashTable;
        std::vector< int > m_hashCollisions;
    };

    class HashMapper
    {
    public:
        HashMapper( const utils::BasePairSequence & paddedHapltypeSequence,
                    const int kmerSize,
                    const int haplotypePadding );

        std::vector< std::size_t > mapSequence( const echidna::utils::BasePairSequence & seq,
                                                boost::optional< std::size_t > hintPosition ) const;

    private:
        KmerMatches m_kmerMatches;
    };
}
}

#endif
