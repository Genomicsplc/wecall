// All content Copyright (C) 2018 Genomics plc
//
// Implementation of a global-local version of the Needleman-Wunsch algorithm using SSE2 instructions.
//
// Gerton Lunter, 20/10/2014
//

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <emmintrin.h>
#include "common.hpp"
#include "align.hpp"
#include "mmHelpers.hpp"

namespace wecall
{
namespace alignment
{

    int needlemanWunschAlignment( std::string::const_iterator haplotype,
                                  std::string::const_iterator readSeq,
                                  const char * readQual,
                                  const unsigned int haplotypeLength,
                                  const unsigned int readLength,
                                  const unsigned short gapextend,
                                  const unsigned short nucprior,
                                  const localGapOpenPenalties_t & localgapopen,
                                  char * aln1,
                                  char * aln2,
                                  int * const firstpos,
                                  const int o )
    {

        /**********************************************************************************************************

        Synopsis:
        =========

        This is an implementation of the Needleman-Wunsch alignment algorithm for a read (readSeq) against a
        piece of genomic sequence (haplotype).  It returns an alignment (if aln1, aln2 != NULL), the minimum score,
        and the position into haplotype that the first (aligning) base of readSeq aligns to.

        The nucprior score penalizes insertions as the alignment score is used as a likelihood of the read
        conditional on the reference, so that any deleted sequence is already conditioned on but any inserted
        sequence is not.

        The algorithm implements a "global-local" version of Needleman-Wunsch, where the full read (readSeq) is
        aligned against an arbitrary subsequence of the genome (haplotype).  Put differently, the initial and
        terminal deletions are not penalized.

        The localgapopen array is of length haplotypeLength, and contains phred scores for opening an insertion or
        deletion at any position in the reference sequence.  A penalty localgapopen[i + o] is charged for deletions
        starting at position i in haplotype, and insertions starting just after base i in haplotype.

        Deletions can immediately follow insertions, but not vice versa -- D->I transitions are forbidden.  This
        ensures that a particular localgapopen penalty can be charged at most once.


        Derivation:
        ===========

        The algorithm can be summarized as

            M(x,y) = min( M(x-1,y-1), I(x-1,y-1), D(x-1,y-1) ) + mismatch(x,y)
            D(x,y) = min( D(x-1,y) + gapextend, M(x-1,y) + gapopen, I(x-1,y) + gapopen )
            I(x,y) = min( I(x,y-1) + gapextend, M(x,y-1) + gapopen ) + nucprior

        where x = 0...haplotypeLength-1 is the genomic position, and y = 0...haplotypeLength-1 is the position in the
        read.  The implementation uses a change of variables,

            s   =   x + y
            d^e =   (x - y)/2
            d^o =   (x - y - 1)/2

        where d^e and d^o are used whenever s is even or odd, respectively.  The M,D,I tables are split into
        two types depending on the parity of  s ,  denoted M^e, D^e, I^e and vice versa, and equations for these
        transformed tables can be derived from the original equations, e.g.

            M^o(s,d) = min( M^o(s-2,d), I^o(s-2,d), D^o(s-2,d) ) + mismatch(s,d)
            D^o(s,d) = min( D^e(s-1,d)   + gapextend, M^e(s-1,d)   + gapopen, I^e(s-1,d) + gapopen )
            I^o(s,d) = min( I^e(s-1,d+1) + gapextend, M^e(s-1,d+1) + gapopen ) + nucprior

        and for M^e

            M^e(s,d) = min( M^e(s-2,d), I^e(s-2,d), D^e(s-2,d) ) + mismatch(s,d)
            D^e(s,d) = min( D^o(s-1,d-1) + gapextend, M^o(s-1,d-1) + gapopen, I^o(s-1,d-1) + gapopen )
            I^e(s,d) = min( I^o(s-1,d)   + gapextend, M^o(s-1,d)   + gapopen ) + nucprior


        Implementation:
        ===============

        1. diagonal banding and SIMD

        The implementation constrains the dynamic programming table to the diagonal  0 <= d^{e,o} < 8 .  This
        allows indels of a maximum size of 14 bp due to the factor 2 in the definition of d^{e,o}.  The main
        reason for constrainng the DP table this way is that the M/D/I arrays are all 8 entries long, allowing
        SIMD instructions to quickly compute the recursion.  Most operations are done by vectorized add and
        minimum instructions; whenever the  d  coordinate changes by +-1, a shift left/right of the vector of
        8 entries is required.  To make sure that all sequence letters can potentially participate in an
        alignment, the longer sequence should be no more than 15 letters longer than the shorter.  This
        implementation requires the longer sequence (haplotype) to be exactly 15 letters longer than the shorter.

        2. in-register DP table

        The algorithm uses the standard low-memory trick of only storing the last few columns (here, diagonals)
        of the DP table.  Because the transformed recursions look back up to 2 positions (in the  s  coordinate)
        the last two diagonals are remembered, requiring 6 vectors.  Since 64-bit processors have 16 XMM registers,
        the entire DP table can be held in registers.  This is the second reason for this implementation's speed.

        3. traceback

        To allow traceback while avoiding storing the entire DP table, a smaller table of pointers is stored
        instead.  Scores in the M, D and I vectors are stored in the upper 14 bits of each 16-bit word, while
        the lower 2 bits store a signature identifying the vector as M, D or I.  The minimum operation carries
        along this signature, so that after one step in the recursion the lower bits denote where the scores
        came from.  The 3 x 2 bits from the 3 vectors are combined into a single vector and stored for traceback.

        4. Initialization

        Initialization is done by loading entry i of the match vector with a score i*offset, pre-loading the
        unused entries in the sequence vector with values guaranteeing a mismatch, and pre-loading unused entries
        in the base quality vector with negative penalties -offset.  The unused entries generate mismatches and
        lower the score, so that precisely when an entry is first used it contains score 0.  This implements the
        penalty-free initial "deletion" of the global-local Needleman-Wunsch algorithm.

        5. Word extraction

        To implement penalty-free deletion at the other end, the match score is extracted when the  y  coordinate
        hits readLength-1.  In the transformed coordinates, this means extracting a single word from a variable position
        within the vector at the last 8 iterations.  The SSE2 instruction set has a word extraction instruction,
        but the index is required to be an immediate.  Instead we use a mask and the MASKMOVDQU instruction that
        stores only those bytes identified by the mask in a memory location, and slide the mask during the last 8
        iterations.

        6. Traceback

        This is entirely standard, after converting back to  x  and  y  coordinates.


        **********************************************************************************************************/

        // bit labels / positions to encode traceback pointers
        constexpr int match_id = 0;
        constexpr int delete_id = 1;
        constexpr int insert_id = 3;

        // check the lengths of the sequences
        assert( haplotypeLength == readLength + 15 );

        // make sure that the initial and final special-case code don't interfere
        assert( haplotypeLength > 8 );

        // check the signs of the penalties
        assert( gapextend > 0 );
        assert( nucprior >= 0 );

        const int backtrace = ( aln1 != nullptr );
        constexpr int init_offset = 2048;
        const short infty = 0x7fff - 4 * gapextend - 4;

        // set the phred score charged for aligning against an N character.  In a probabilistic diploid this
        // is the probability of drawing a particular nucleotide from the background distribution.  This is
        // the same probability as for inserting a particular nucleotide, so use that phred score.
        const int nscore = 4 * nucprior;

        const short_array8 gapextend_vec( short( 4 * gapextend ) );
        const short_array8 nucprior_vec( short( 4 * nucprior ) );
        const short_array8 three_vec( short( 3 ) );

        auto localgapopenStart = localgapopen.begin() + o;

        // this (stack-allocated) variable is used to extract words from 128-bit arrays
        short_array8 wordvector;

        // the pointers used for the backtrace algorithm
        // (For the traceback algorithm, the values in traceback_ptrs are accessed through _wordbackpointers;
        //  according to the GCC documentation, adding __may_alias__ to the latter's definition should ensure
        //  this works, however it proved to be necessary to also give traceback_ptrs this attribute.  Bit worrying...)
        std::vector< short_array8 > traceback_ptrs( 2 * ( haplotypeLength + 8 ) );
        int16_t * word_traceback_ptrs = (int16_t *)&traceback_ptrs[0];

        // initialization, to -32768 (=(int16_t)0x8000) which encodes score 0, plus offsets that are annulled by
        // mismatches
        constexpr int x8000 = 0x8000;
        short_array8 mat_even( [init_offset]( int i )
                               {
                                   return short( 0x8000 + i * init_offset );
                               } );
        auto mat_odd = mat_even;
        short_array8 ins_even( infty );
        short_array8 del_even( infty );
        short_array8 ins_odd( infty );
        short_array8 del_odd( infty );

        // initialize seq1vector with the first 8 characters of sequence 1
        short_array8 seq1vector( [haplotype]( int i )
                                 {
                                     return haplotype[i];
                                 } );
        short_array8 qual1vector( [haplotype, nscore, infty]( int i )
                                  {
                                      return haplotype[i] == constants::gapChar ? nscore : infty;
                                  } );
        short_array8 seq2vector = ins_even;                 // this ensures mismatches for out-of-bound chars
        short_array8 qual2vector( short( -init_offset ) );  // this ensures cells are initialized correctly
        short_array8 gapopen_vec( [localgapopenStart, o]( int i )
                                  {
                                      return 4 * localgapopenStart[i];
                                  } );

        //
        // main loop
        //

        unsigned int s;
        int minscore = 99999;
        int s_backtrace = -1;

        for ( s = 0; s < 2 * readLength + 14; s += 2 )
        {

            /*
            printxmm("gapopen   ",gapopen_vec);
            printxmm("qual1vec  ",qual1vector);
            printf("s=%i ",s-2);printxmm("besteven",mat_even);
            printf("s=%i ",s-1);printxmm("mat_odd ",mat_odd);
            printf("s=%i ",s-1);printxmm("ins_odd ",ins_odd);
            printf("s=%i ",s-1);printxmm("del_odd ",del_odd);
            */

            // haplotype is current; readSeq needs updating
            seq2vector = shift_left( seq2vector );
            qual2vector = shift_left( qual2vector );

            if ( s / 2 < readLength )
            {
                seq2vector[0] = readSeq[s / 2];
                qual2vector[0] = 4 * readQual[s / 2];
            }

            //
            // first handle the case of s even
            //

            // extract minimum score in mat_even; the calculation of the minimum across mat/del/ins is done at the end
            // of the current block to free up some registers; and it is not necessary at the very first iteration.
            if ( s / 2 >= readLength )
            {
                if ( mat_even[s / 2 - readLength] < minscore )
                {
                    minscore = mat_even[s / 2 - readLength];
                    s_backtrace = s;
                }
            }

            // the actual recursion
            mat_even = mat_even + min( qual1vector, andnot( cmpeq( seq2vector, seq1vector ), qual2vector ) );

            ins_even = min( mat_odd + gapopen_vec, ins_odd + gapextend_vec ) + nucprior_vec;

            del_even = min( min( mat_odd, ins_odd ) + shift_right( gapopen_vec ), del_odd + gapextend_vec );

            del_even = shift_left( del_even );
            del_even[0] = infty;

            // doing this calculation here frees up the registers storing ins_odd and del_odd
            mat_odd = min( mat_odd, min( ins_odd, del_odd ) );

            /*
            printf("s=%i ",s-1); printxmm("bestodd ",mat_odd);
            printf("s=%i ",s);printxmm("mat_even ",mat_even);
            printf("s=%i ",s);printxmm("ins_even ",ins_even);
            printf("s=%i ",s);printxmm("del_even ",del_even);
            */

            // get back-pointers and store them.
            //
            //    traceback_ptrs[ s ] =
            //           (three_vec && mat_even) || ((three_vec && ins_even) << 2*insert_id) || ((three_vec && del_even)
            //           << 2*delete_id)
            //
            // where && || and << denote the 128-bit versions of these operations.
            if ( backtrace )
            {
                ins_odd = three_vec & ins_even;
                ins_odd = ins_odd << 2 * insert_id;
                del_odd = three_vec & mat_even;
                ins_odd = ins_odd | del_odd;
                del_odd = three_vec & del_even;
                del_odd = del_odd << 2 * delete_id;
                traceback_ptrs[s] = ins_odd | del_odd;

                // set state labels
                mat_even = andnot( three_vec, mat_even );
                del_even = andnot( three_vec, del_even ) | ( three_vec >> 1 );
                ins_even = andnot( three_vec, ins_even ) | three_vec;
            }

            //
            // now handle case of s odd (but don't bother actually changing the s variable)
            //

            // update haplotype
            if ( s / 2 + 8 < haplotypeLength )
            {
                ins_odd = shift_right( seq1vector );
                seq1vector = ins_odd;
                seq1vector[7] = haplotype[s / 2 + 8];
                del_odd = shift_right( qual1vector );
                qual1vector = del_odd;
                qual1vector[7] = haplotype[s / 2 + 8] == constants::gapChar ? nscore : infty;
                ins_odd = shift_right( gapopen_vec );
                gapopen_vec = ins_odd;
                gapopen_vec[7] = 4 * localgapopenStart[s / 2 + 8];
            }
            else
            {
                seq1vector = shift_right( seq1vector );
                qual1vector = shift_right( qual1vector );
                gapopen_vec = shift_right( gapopen_vec );
            }

            // at this point, extract minimum score.  Referred-to position must be y==readLength-1, so that current
            // position has y==readLength; i==0 so d=0 and y=s/2
            if ( s / 2 >= readLength )
            {
                if ( mat_odd[s / 2 - readLength] < minscore )
                {
                    minscore = mat_odd[s / 2 - readLength];
                    s_backtrace = s + 1;
                }
            }

            // the actual recursion
            mat_odd = mat_odd + min( qual1vector, andnot( cmpeq( seq2vector, seq1vector ), qual2vector ) );

            del_odd = min( min( mat_even, ins_even ) + gapopen_vec, del_even + gapextend_vec );

            ins_odd = min( ( shift_right( mat_even ) + gapopen_vec ), shift_right( ins_even ) + gapextend_vec ) +
                      nucprior_vec;

            ins_odd[7] = infty;

            // doing this calculation here frees up the registers storing ins_even and del_even
            mat_even = min( mat_even, min( ins_even, del_even ) );

            // get back-pointers by calculating
            //
            //    traceback_ptrs[ s+1 ] = (three_vec && mat_odd) || ((three_vec && ins_odd) << 2*insert) || ((three_vec
            //    && del_odd) << 2*delete)
            //
            if ( backtrace )
            {
                ins_even = three_vec & ins_odd;
                ins_even = ins_even << 2 * insert_id;
                del_even = three_vec & mat_odd;
                ins_even = ins_even | del_even;
                del_even = three_vec & del_odd;
                del_even = del_even << 2 * delete_id;
                traceback_ptrs[s + 1] = ins_even | del_even;

                // set state labels
                mat_odd = andnot( three_vec, mat_odd );
                del_odd = andnot( three_vec, del_odd ) | ( three_vec >> 1 );
                ins_odd = andnot( three_vec, ins_odd ) | three_vec;
            }
        }

        // do last iteration -- since only the mat_* values are used, and mat_even is already calculated,
        // this boils down to taking minima across mat/del/ins to calculate mat_odd
        mat_odd = min( mat_odd, min( ins_odd, del_odd ) );

        // put down last scores
        if ( mat_even[s / 2 - readLength] < minscore )
        {
            minscore = mat_even[s / 2 - readLength];
            s_backtrace = s;
        }
        if ( mat_odd[s / 2 - readLength] < minscore )
        {
            minscore = mat_odd[s / 2 - readLength];
            s_backtrace = s + 1;
        }

        if ( backtrace )
        {
            traceback_ptrs[s] = mat_even;
            traceback_ptrs[s + 1] = mat_odd;
        }

        // move to 0-based integers
        minscore += x8000;

        // if no backtrace required, we're done
        if ( not backtrace )
        {
            return minscore >> 2;
        }

        // at this point, s_backtrace points to a dummy Match state representing the lowest score
        int d = s_backtrace / 2 - readLength;
        unsigned int state = ( word_traceback_ptrs[8 * s_backtrace + d] >> ( 2 * match_id ) ) & 3;
        int i = 0;
        int x = s_backtrace - readLength;
        int y = readLength;

        // move to the first actual alignment state
        s_backtrace -= 2;

        // build the alignment backwards-to-front
        while ( y > 0 )
        {
            const unsigned int previous_state = ( word_traceback_ptrs[8 * s_backtrace + d] >> ( 2 * state ) ) & 3;
            if ( state == match_id )
            {
                aln1[i] = haplotype[--x];
                aln2[i] = readSeq[--y];
                s_backtrace -= 2;
            }
            else if ( state == insert_id )
            {
                d += s_backtrace & 1;
                aln1[i] = '-';
                aln2[i] = readSeq[--y];
                s_backtrace -= 1;
            }
            else
            {
                d -= ( s_backtrace - 1 ) & 1;
                aln1[i] = haplotype[--x];
                aln2[i] = '-';
                s_backtrace -= 1;
            }
            state = previous_state;
            i++;
        }
        aln1[i] = aln2[i] = 0;

        if ( firstpos != nullptr )
        {
            *firstpos = x;
        }

        // reverse the alignment, by swapping pairs of columns and moving inwards
        for ( auto it = aln1, rit = aln1 + i - 1; it < rit; ++it, --rit )
        {
            std::swap( *it, *rit );
        }
        for ( auto it = aln2, rit = aln2 + i - 1; it < rit; ++it, --rit )
        {
            std::swap( *it, *rit );
        }

        // and done
        return minscore >> 2;
    }
}
}
