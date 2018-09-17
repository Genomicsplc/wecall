// All content Copyright (C) 2018 Genomics plc
#ifndef ERROR_CORRECTION_PARAMETERS_HPP
#define ERROR_CORRECTION_PARAMETERS_HPP

#include <vector>
#include <utility>

#include "common.hpp"

namespace echidna
{
namespace corrector
{
    struct ErrorCorrectionParameters
    {
        using weightedDistance = std::pair< int, double >;  ///< slippage distance and probability that it occurs

        const double pm = 0.8;       ///< probability of nucleotide match in error state
        const double ps = 0.8;       ///< proportion of nucleotide mismatches explained by slippage
        const double ptrue = 0.5;    ///< proportion of kmers in error state that have proper q scores
        const double perr = 0.0005;  ///< (prior) per-nuc probability of read turning into error state
        const double pref = 0.95;    ///< probability of true kmer being reference
        const std::vector< weightedDistance > distances = {{-1, 0.5}, {1, 0.5}};

        const char qualityFloor = 0x05;  ///<Avoid accumulating weak and systematic evidence.
        const char floorTo =
            constants::minAllowedQualityScore;  ///< Gives nearly equal biases amongst the 4 base pairs.

        // TODO(semen): check that probs in distances add up to unity
        // TODO(semen): check that max abs distances == padding
    };
}  // namespace corrector
}  // namespace echidna

#endif
