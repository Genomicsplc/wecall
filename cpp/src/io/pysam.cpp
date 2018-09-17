// All content Copyright (C) 2018 Genomics plc
#include <ctype.h>
#include <assert.h>

// Declare the samtools headers as extern
extern "C" {
#include "samtools/sam.h"
#include "samtools/khash.h"
#include "samtools/ksort.h"
}

#include "io/pysam.hpp"

// #######################################################
// utility routines to avoid using callbacks in bam_fetch
// taken from bam_index.c
// The order of the following declarations is important.
// #######################################################

struct pair64_t
{
    uint64_t u, v;
};

#define pair64_lt( a, b ) ( ( a ).u < ( b ).u )

struct bam_binlist_t
{
    uint32_t m, n;
    pair64_t * list;
};

struct bam_lidx_t
{
    int32_t n, m;
    uint64_t * offset;
};

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"

// cppcheck-suppress variableScope
// cppcheck-suppress cstyleCast
KSORT_INIT( my_off, pair64_t, pair64_lt )
// cppcheck-suppress cstyleCast
KHASH_MAP_INIT_INT( my_i, bam_binlist_t )

#pragma clang diagnostic pop

struct __bam_index_t
{
    int32_t n;
    khash_t( my_i ) * *index;
};

// standin for bam_destroy1 in bam.h
// deletes all variable length data
void pysam_bam_destroy1( bam1_t * b )
{
    if ( b == NULL )
    {
        return;
    }
    if ( b->data != NULL )
    {
        free( b->data );
    }
    free( b );
}

//-------------------------------------------------------------------------------------------------

uint32_t * pysam_bam1_cigar( const bam1_t * b ) { return (uint32_t *)( b->data + b->core.l_qname ); }

//-------------------------------------------------------------------------------------------------

uint8_t * pysam_bam1_seq( const bam1_t * b ) { return (uint8_t *)( b->data + b->core.n_cigar * 4 + b->core.l_qname ); }

//-------------------------------------------------------------------------------------------------

uint8_t * pysam_bam1_qual( const bam1_t * b )
{
    return (uint8_t *)( b->data + b->core.n_cigar * 4 + b->core.l_qname + ( b->core.l_qseq + 1 ) / 2 );
}

//-------------------------------------------------------------------------------------------------
//
//
//
// Iterator implementation follows
//
//
// functions defined in bam_index.c

extern "C" pair64_t * get_chunk_coordinates( const bam_index_t * idx, int tid, int beg, int end, int * cnt_off );

//-------------------------------------------------------------------------------------------------

static inline int is_overlap( uint32_t beg, uint32_t end, const bam1_t * b )
{
    uint32_t rbeg = b->core.pos;
    uint32_t rend = b->core.n_cigar ? bam_calend( &b->core, bam1_cigar( b ) ) : b->core.pos + 1;
    return ( rend > beg && rbeg < end );
}

//-------------------------------------------------------------------------------------------------

struct __bam_fetch_iterator_t
{
    bam1_t * b;
    pair64_t * off;
    int n_off;
    uint64_t curr_off;
    int curr_chunk;
    bamFile fp;
    int tid;
    int beg;
    int end;
    int n_seeks;
};

//-------------------------------------------------------------------------------------------------

bam_fetch_iterator_t * bam_init_fetch_iterator( bamFile fp, const bam_index_t * idx, int tid, int beg, int end )
{
    // iterator contains current alignment position
    //      and will contain actual alignment during iterations
    bam_fetch_iterator_t * iter = (bam_fetch_iterator_t *)calloc( 1, sizeof( bam_fetch_iterator_t ) );
    iter->b = (bam1_t *)calloc( 1, sizeof( bam1_t ) );

    // list of chunks containing our alignments
    iter->off = get_chunk_coordinates( idx, tid, beg, end, &iter->n_off );

    // initialise other state variables in iterator
    iter->fp = fp;
    iter->curr_chunk = -1;
    iter->curr_off = 0;
    iter->n_seeks = 0;
    iter->tid = tid;
    iter->beg = beg;
    iter->end = end;
    return iter;
}

//-------------------------------------------------------------------------------------------------

bam1_t * bam_fetch_iterate( bam_fetch_iterator_t * iter )
{
    if ( not iter->off )
    {
        return 0;
    }

    // iterate through all alignments in chunks
    for ( ;; )
    {
        // then jump to the next chunk
        if ( iter->curr_off == 0 || iter->curr_off >= iter->off[iter->curr_chunk].v )
        {
            if ( iter->curr_chunk == iter->n_off - 1 )
            {
                // no more chunks
                break;
            }
            if ( iter->curr_chunk >= 0 )
            {
                // otherwise bug
                assert( iter->curr_off == iter->off[iter->curr_chunk].v );
            }
            if ( iter->curr_chunk < 0 || iter->off[iter->curr_chunk].v != iter->off[iter->curr_chunk + 1].u )
            {
                // not adjacent chunks; then seek
                bam_seek( iter->fp, iter->off[iter->curr_chunk + 1].u, SEEK_SET );
                iter->curr_off = bam_tell( iter->fp );
                ++iter->n_seeks;
            }
            ++iter->curr_chunk;
        }
        if ( bam_read1( iter->fp, iter->b ) > 0 )
        {
            iter->curr_off = bam_tell( iter->fp );
            if ( iter->b->core.tid != iter->tid || iter->b->core.pos >= iter->end )
            {
                // no need to proceed
                break;
            }
            else if ( is_overlap( iter->beg, iter->end, iter->b ) )
            {
                // func(iter->b, data);
                return iter->b;
            }
        }
        else
        {
            return 0;  // end of file
        }
    }
    return 0;
}

//-------------------------------------------------------------------------------------------------

void bam_cleanup_fetch_iterator( bam_fetch_iterator_t * iter )
{
    bam_destroy1( iter->b );
    free( iter->off );
}

//-------------------------------------------------------------------------------------------------
