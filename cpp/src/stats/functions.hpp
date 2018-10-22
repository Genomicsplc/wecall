// All content Copyright (C) 2018 Genomics plc
#ifndef STATS_FUNCTIONS_HPP
#define STATS_FUNCTIONS_HPP

/// @file
/// @author Andrew Rimmer
/// @version 0.0
///
/// @section LICENSE
/// I have no idea what kind of license we are using for this..
/// but we should put the license here
///
/// @section Description
/// Statistical functions for use in weCall.

#include "common.hpp"

namespace wecall
{
namespace stats
{
    /// Compute the cumulative distribution function for the beta-binomial
    /// distribution, which is a Binomial with a Beta prior for the special case used in reference calls
    /// (alpha == beta, k=0)
    ///
    /// @param n The number of trials
    /// @param alpha Parameter for the Beta distribution used as prior
    /// @return The cumulative value of the distribution for k=0
    double betaBinomialCDFForReferenceCalls( int n, const double alpha );

    /// Compute the cumulative distribution function for the beta-binomial
    /// distribution, which is a Binomial with a Beta prior.
    ///
    /// @param k The number of successes
    /// @param n The number of trials
    /// @param alpha Parameter for the Beta distribution used as prior
    /// @param beta Parameter for the Beta distribution used as prior
    /// @return The cumulative value of the distribution up to k
    double betaBinomialCDFStable( const int k, const int n, const int alpha, const int beta );
    double betaBinomialCDF( const int k, const int n, const double alpha, const double beta );

    /// The natural logrithm of the Beta function i.e. ln[ Gamma(x) * Gamma(y) / Gamma(x + y) ]
    ///
    /// @param x
    /// @param y
    /// @return ln(beta(x,y))
    double logBetaFunction( const double x, const double y );

    /// Calculate the value of the generalised hypergeometric function 3F2. This is
    /// currently only used for computing the Beta-binomial CDF, which is why the
    /// parameters are n,k,alpha,beta.
    ///
    /// @param n Number of trials
    /// @param k Number of successes
    /// @param alpha parameter of beta distribution
    /// @param beta  parameter of beta distribution
    /// @return value of 3F2
    double threeFtwo( const double a2, const double a3, const double b1, const double b2 );

    constexpr double phredCoefficient = -10.0;
    /// Convert probability into Phred score.
    ///
    /// @param probability Input probability
    /// @return Output Phred score
    phred_t toPhredQ( const double probability );

    /// Convert Phred score into probability.
    ///
    /// @param phredQ Input Phred score
    /// @return Output probability
    double fromPhredQ( const phred_t phredQ );

    unsigned int roundPhred( const phred_t phred );

    double variantSupportPerRead( const double prior, const double posterior, const int64_t variantSupportCount );
}
}

#endif
