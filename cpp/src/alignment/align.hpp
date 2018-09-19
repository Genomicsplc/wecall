// All content Copyright (C) 2018 Genomics plc
#ifndef ALIGN_HPP
#define ALIGN_HPP

#include <vector>
#include <emmintrin.h>
#include <stdio.h>
#include <assert.h>

namespace echidna
{
namespace alignment
{
    using errorModel_t = std::vector< int16_t >;
    using localGapOpenPenalties_t = errorModel_t;

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
                                  const int o  ///< offset
                                  );
}
}

#endif
