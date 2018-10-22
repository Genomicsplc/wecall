// All content Copyright (C) 2018 Genomics plc
#include "mapping/hashMapper.hpp"
#include "utils/exceptions.hpp"
#include "utils/logging.hpp"
#include "alignment/galign.hpp"
#include "utils/bestScoreSelector.hpp"

#include <stdexcept>
#include <cmath>
#include <iostream>
#include <map>
#include <cassert>

namespace wecall
{
namespace mapping
{

    HashFunction::HashFunction( const utils::BasePairSequence & sequence, unsigned int kmerSize )
        : m_bitShift( kmerSize * 2 - 3 )
    {
        WECALL_ASSERT( kmerSize > 1, "Cannot operate with kmers of size 1" );
        WECALL_ASSERT( sequence.size() >= kmerSize, "Sequence is too short for hashing. Seq = " + sequence.str() );

        const auto maxSeqLen = 0x1u << 2 * kmerSize;
        assert( maxSeqLen == std::pow( 4, kmerSize ) );
        WECALL_ASSERT( sequence.size() < maxSeqLen, "Sequence is too long for hashing. Length = " +
                                                         std::to_string( sequence.size() ) + " and max allowed = " +
                                                         std::to_string( maxSeqLen ) );

        m_currentHash = 0;
        auto begin = sequence.cbegin();
        auto end = begin + kmerSize - 1;
        int shift = 1;

        for ( auto it = begin; it != end; ++it, shift += m_nBits )
        {
            const int hashAdd = ( ( *it ) & m_bitMask );
            m_currentHash |= ( hashAdd << shift );
        }
    }

    int HashFunction::next( const char newBase )
    {
        const int hashAdd = newBase & m_bitMask;
        m_currentHash = ( m_currentHash >> m_nBits ) | ( hashAdd << m_bitShift );
        return m_currentHash;
    }

    KmerMatches::KmerMatches( const utils::BasePairSequence & paddedHaplotypeSequence,
                              const int kmerSize,
                              const int haplotypePadding )
        : m_paddedHaplotypeSequenceLength( static_cast< int >( paddedHaplotypeSequence.size() ) ),
          m_kmerSize( kmerSize ),
          m_haplotypePadding( haplotypePadding ),
          m_hashSize( std::pow( 4, m_kmerSize ) ),
          m_hashTable( m_hashSize, constants::endOfChain ),
          m_hashCollisions( paddedHaplotypeSequence.size(), constants::endOfChain )
    {
        this->indexSequence( paddedHaplotypeSequence );
    }

    std::vector< std::size_t > KmerMatches::countKmerMatches(
        const wecall::utils::BasePairSequence & readSequence ) const
    {
        WECALL_ASSERT( readSequence.size() <= m_paddedHaplotypeSequenceLength,
                        "Require padded Haplotype sequence to be longer than any chosen read-sequence" );

        const int readSeqLen = static_cast< int >( readSequence.size() );

        HashFunction hashFunction( readSequence, m_kmerSize );

        const auto allowablePositions = this->allowableStartPositionsInHaplotype( readSeqLen );

        std::vector< std::size_t > counts( int64_to_sizet( allowablePositions.end() ), 0 );

        auto readSequenceIt = readSequence.cbegin() + m_kmerSize - 1;
        for ( std::size_t index = 0, index_end = readSeqLen - ( m_kmerSize - 1 ); index < index_end;
              ++index, ++readSequenceIt )
        {
            const auto base = *readSequenceIt;

            const int hashVal = hashFunction.next( base );

            auto hapIdx = m_hashTable[hashVal];

            while ( this->isHashHit( hapIdx ) )
            {
                const auto pos = hapIdx - index;

                if ( allowablePositions.contains( pos ) )
                {
                    ++counts[pos];
                }

                hapIdx = m_hashCollisions[hapIdx];
            }
        }

        return counts;
    }

    void KmerMatches::indexSequence( const utils::BasePairSequence & haplotypeSequence )
    {
        HashFunction hashFunction( haplotypeSequence, m_kmerSize );

        for ( std::size_t index = 0; index < haplotypeSequence.size() - ( m_kmerSize - 1 ); ++index )
        {
            const auto base = haplotypeSequence[index + m_kmerSize - 1];
            int hashVal = hashFunction.next( base );

            if ( m_hashTable[hashVal] == constants::endOfChain )
            {
                // No entry in the hash table yet with this hashVal
                m_hashTable[hashVal] = index;
            }
            else
            {
                // Collision
                int j = m_hashTable[hashVal];
                int count = 2;  // the existing entry, plus the new one

                while ( j != constants::repeatChain && m_hashCollisions[j] != constants::endOfChain )
                {
                    j = m_hashCollisions[j];
                    ++count;
                }

                if ( count > constants::maxRepeatCount )
                {
                    // remove the entire chain
                    m_hashTable[hashVal] = constants::repeatChain;
                }
                else
                {
                    // add new entry, unless this is a repeat chain
                    if ( j != constants::repeatChain )
                    {
                        m_hashCollisions[j] = index;
                    }
                }
            }
        }
    }

    utils::Interval KmerMatches::allowableStartPositionsInHaplotype( const int64_t readSeqLen ) const
    {
        return alignment::allowableStartPositionsForAlignment( m_paddedHaplotypeSequenceLength, readSeqLen,
                                                               m_haplotypePadding );
    }

    HashMapper::HashMapper( const utils::BasePairSequence & paddedHapltypeSequence,
                            const int kmerSize,
                            const int haplotypePadding )
        : m_kmerMatches( paddedHapltypeSequence, kmerSize, haplotypePadding )
    {
    }

    //-----------------------------------------------------------------------------------------

    std::vector< std::size_t > HashMapper::mapSequence( const utils::BasePairSequence & sequence,
                                                        boost::optional< std::size_t > hintPosition ) const
    {
        // Identify and collect the top bins in the histogram
        const auto counts = m_kmerMatches.countKmerMatches( sequence );
        const double maxHighestLowestRatio = 5.0;
        std::vector< std::size_t > matches = utils::indiciesWithHighestValues(
            counts, counts.size(), constants::fractionOfKmersToConsider, maxHighestLowestRatio );

        if ( matches.empty() and hintPosition )
        {
            const auto allowablePositions = m_kmerMatches.allowableStartPositionsInHaplotype( sequence.size() );
            std::size_t allowableStart = int64_to_sizet( allowablePositions.start() );

            if ( allowablePositions.contains( hintPosition.get() ) )
            {
                matches.push_back( hintPosition.get() );
            }
            else if ( hintPosition < allowableStart )
            {
                matches.push_back( allowableStart );
            }
            else
            {
                matches.push_back( int64_to_sizet( allowablePositions.end() - 1 ) );
            }
        }

        return matches;
    }

    //-----------------------------------------------------------------------------------------
}
}
