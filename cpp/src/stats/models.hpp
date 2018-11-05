// All content Copyright (C) 2018 Genomics plc
#ifndef STATS_MODELS_HPP
#define STATS_MODELS_HPP

/// @file
/// @author Andrew Rimmer
/// @version 0.0
///
/// @section LICENSE
/// I have no idea what kind of license we are using for this..
/// but we should put the license here
///
/// @section Description
/// Utility functions related to some of the probabilistic models used
/// in weCall
#include <vector>
#include <algorithm>
#include <numeric>

#include "common.hpp"

namespace wecall
{
namespace stats
{
    /// Use a beta-binomial distribution to diploid the probability that all reads
    /// at a particular locus are from the same haplotype, i.e. the same copy of a chromosome.
    /// Distribution parameters alpha and beta have been tuned empirically, based on
    /// whole genome data from Illumina.
    ///
    /// @param nReads The number of reads covering the locus
    /// @return The probablity that all reads originate from the same haplotype
    double probAllReadsFromSingleCopyOfDiploid( int nReads );

    /// Use a beta-binomial distribution to diploid the probability that the number of reads that support
    /// a given variant are sufficient to support a heterozygous call of that variant (allele bias P-value).
    /// Distribution parameters alpha and beta have been tuned empirically, based on
    /// whole genome data from Illumina.
    ///
    /// @param totalCoverage The total number of reads covering the variant locus.
    /// @param varCoverage The number of those reads that support the variant.
    /// @return The allele bias P-value.
    double probSufficientVarCoverageToSupportHet( int totalCoverage, int varCoverage );

    /// Use a beta-binomial distribution to diploid the probability that support for a variant
    /// across strands is consistent with coverage across strands (strand bias P-value).
    /// Distribution parameters alpha and beta have been tuned empirically, based on
    /// whole genome data from Illumina.
    ///
    /// @param totalFwdCoverage The total number of forward reads covering the variant locus.
    /// @param varFwdCoverage The number of forward reads that support the variant.
    /// @param totalRevCoverage The total number of reverse reads covering the variant locus.
    /// @param varRevCoverage The number of reverse reads that support the variant.
    /// @return The strand bias P-value.
    double probVarSupportNotBiasedByStrand( int totalFwdCoverage,
                                            int varFwdCoverage,
                                            int totalRevCoverage,
                                            int varRevCoverage );

    template < typename Type >
    double rootMeanSquare( const std::vector< Type > & values )
    {
        if ( values.size() == 0 )
        {
            return std::numeric_limits< double >::quiet_NaN();
        }
        else
        {
            auto square = []( double sum, Type a )
            {
                const double double_a = static_cast< double >( a );
                return sum + double_a * double_a;
            };
            return sqrt( std::accumulate( values.cbegin(), values.cend(), 0.0, square ) /
                         static_cast< double >( values.size() ) );
        }
    }
}
}

#endif
