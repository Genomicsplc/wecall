// All content Copyright (C) 2018 Genomics plc
#ifndef HAPLOTYPE_LIKELIHOODS_HPP
#define HAPLOTYPE_LIKELIHOODS_HPP

#include "utils/matrix.hpp"
#include "variant/haplotype.hpp"
#include "io/readRange.hpp"
#include "alignment/galign.hpp"

namespace wecall
{
namespace caller
{
    namespace errorModels
    {
        const wecall::alignment::errorModel_t illuminaErrorModel = {
            45, 42, 41, 39, 37, 32, 28, 23, 20, 19, 17, 16, 15, 14, 13, 12, 11, 11, 10, 9, 9, 8, 8, 7, 7,
            7,  6,  6,  6,  5,  5,  5,  4,  4,  4,  3,  3,  3,  3,  2,  2,  2,  2,  2,  1, 1, 1, 1, 1};
    }

    utils::matrix_t computeHaplotypeLikelihoods( const variant::HaplotypeVector & haplotypes,
                                                 const io::RegionsReads & readRange );

    std::vector< double > computeHaplotypeFrequencies( const utils::matrix_t & haplotypeLikelihoods );
}
}

#endif
