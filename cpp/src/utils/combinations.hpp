// All content Copyright (C) 2018 Genomics plc
#ifndef COMBINATIONS_HPP
#define COMBINATIONS_HPP

#include <vector>
#include <set>
#include <climits>

namespace wecall
{
namespace utils
{
    //-------------------------------------------------------------------------------------------------

    /// Return the next permutation of N 1 bits, where 'next' means next in lexicographical order.
    /// For example, if N is 3 and the bit pattern is 00010011, the next patterns would be
    //  00010101, 00010110, 00011001,00011010, 00011100, 00100011,
    /// and so forth. This is based on an code from the following url:
    /// http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
    ///
    /// @param input An integer, treated as a 32 bit bitset
    /// @return An integer representing the next permutation of the bits in 'input'
    inline std::size_t getNextBitPermutation( const int64_t input )
    {
        const std::size_t temp = ( input | ( input - 1 ) ) + 1;
        const std::size_t output = temp | ( ( ( ( temp & -temp ) / ( input & -input ) ) >> 1 ) - 1 );
        return output;
    }

    //-------------------------------------------------------------------------------------------------

    /// Generates all N-tuples of input. N is the number of objects in each output subset. This method uses
    /// some bit-twiddling to generate all ordered bit permutations where N bits are set in a standard unsigned
    /// integer.
    ///
    /// @param objects The set of objects from which we want to return subsets
    /// @param N The size of each returned subset
    /// @return A vector of subsets of the input objects
    template < typename C, typename T >
    std::vector< C > combinations( const std::vector< T > & objects, const std::size_t N )
    {
        const std::size_t nObj = objects.size();

        assert( nObj <= sizeof( int64_t ) * CHAR_BIT );  // Make sure we have enough bits to play with
        assert( N <= nObj );                             // Don't try to take more objects than there are in the vector

        std::vector< C > output;

        std::size_t first = 0;
        std::size_t last = 0;

        for ( std::size_t i = 0; i < N; ++i )
        {
            first |= ( 1 << i );                // Set first N bits for starting value
            last |= ( 1 << ( nObj - 1 - i ) );  // Set last N bits for final value
        }

        for ( std::size_t bits = first; bits <= last; bits = getNextBitPermutation( bits ) )
        {
            C newComb;

            for ( std::size_t testBit = 1, pos = 0; testBit < ( 1 << nObj ); testBit <<= 1, ++pos )
            {
                if ( testBit & bits )  /// Is this bit set?
                {
                    newComb.insert( newComb.end(), objects[pos] );
                }
            }

            output.push_back( newComb );
        }

        return output;
    }

    //-------------------------------------------------------------------------------------------------
}
}

#endif
