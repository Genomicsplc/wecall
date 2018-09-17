// All content Copyright (C) 2018 Genomics plc
#include "stats/functions.hpp"

#include <cassert>
#include <algorithm>
#include <cmath>
#include "utils/logging.hpp"

namespace echidna
{
namespace stats
{
    // Calculate cumulative probability for a beta-binomial distribution for special case k=0, alpha=beta
    double betaBinomialCDFForReferenceCalls( int n, const double alpha )
    {
        // case: n == k
        if ( n == 0 )
        {
            return 1.0;
        }

        double one_minus_cdf = alpha * n / ( 2 * alpha + n - 1 );
        for ( int i = 0; i < alpha; ++i )
        {
            one_minus_cdf *= ( alpha + i ) / ( alpha + n - 1 + i );
        }

        // bound result by machine precision (kicks in for n=126)
        if ( std::abs( one_minus_cdf ) < 4.9e-14 )
        {
            return 1e-15;
        }
        const double cdf = 1 - one_minus_cdf * threeFtwo( alpha + 1, -n + 1, 2, -alpha - n + 2 );
        return cdf;
    }

    // Calculate cumulative probability for a beta-binomial distribution
    double betaBinomialCDFStable( const int k, const int n, const int alpha, const int beta )
    {
        // Cumulative frequency is guaranteed to be 1.0. 3F2 doesn't calculate for these values.
        if ( k == n )
        {
            return 1.0;
        }

        double threeFtwo_term = threeFtwo( alpha + k + 1, -n + k + 1, k + 2, -beta - n + k + 2 );

        double term1 = 1;
        for ( int i = n - k; i <= n; ++i )
        {
            term1 *= i;
        }

        double term2 = 1;
        for ( int i = alpha; i <= alpha + k; ++i )
        {
            term2 *= i;
        }

        double term3 = 1;
        for ( int i = 1; i <= k + 1; ++i )
        {
            term3 *= i;
        }

        long double term4 = 1.0;
        for ( int i = 0; i <= n - k - 2; ++i )
        {
            term4 *= static_cast< long double >( beta + i ) / static_cast< long double >( alpha + beta + i );
        }
        for ( int i = n - k - 1; i <= n - 1; ++i )
        {
            term4 /= ( alpha + beta + i );
        }

        double one_minus_cdf = term1 * term2 / term3 * term4;
        // return default value once we run into machine precision problems (happens for ref calls at n=126)
        if ( one_minus_cdf < 4.8e-14 )
        {
            return 1e-15;
        }

        double cdf = 1.0 - one_minus_cdf * threeFtwo_term;
        return cdf;
    }

    // Calculate cumulative probability for a beta-binomial distribution
    double betaBinomialCDF( const int k, const int n, const double alpha, const double beta )
    {
        if ( k == n )
        {
            // Cumulative frequency is guaranteed to be 1.0. 3F2 doesn't calculate for these values.
            return 1.0;
        }
        const double numerator = logBetaFunction( beta + n - k - 1, alpha + k + 1 ) +
                                 log( threeFtwo( alpha + k + 1, -n + k + 1, k + 2, -beta - n + k + 2 ) );
        const double denominator = logBetaFunction( alpha, beta ) + logBetaFunction( n - k, k + 2 ) + log( n + 1 );

        const double relativeDiff = ( numerator - denominator ) / ( numerator + denominator );
        double betaBinomialCDF;

        // use exact formula or floor distribution when we reach machine precision
        if ( std::abs( relativeDiff ) > 1e-15 )
        {
            betaBinomialCDF = 1.0 - exp( numerator - denominator );
        }
        else
        {
            betaBinomialCDF = -1e-15 * ( numerator + denominator );
        }

        return betaBinomialCDF;
    }

    // Calculate the natural logarithm of the beta function
    // log[ Gamma(x) * Gamma(y) / Gamma(x + y) ]
    double logBetaFunction( const double x, const double y ) { return lgamma( x ) + lgamma( y ) - lgamma( x + y ); }

    // Calculate the generalised hypergeometric function 3F2 with a1 = 1 and z = 1.
    double threeFtwo( const double a2, const double a3, const double b1, const double b2 )
    {
        double lastTerm = 1.0;
        double sum = 1.0;                               // term for i = 0
        for ( int i = 1; i < std::abs( a3 ) + 1; ++i )  // Limit on i here is arbitrary.  Copied from the old code..
        {
            const double newTerm =
                lastTerm * ( a2 + i - 1.0 ) * ( a3 + i - 1.0 ) / ( ( b1 + i - 1.0 ) * ( b2 + i - 1.0 ) );
            sum += newTerm;
            lastTerm = newTerm;
        }
        return sum;
    }

    // Converts an error probability to a "PHRED" quality score (1.0 -> 0, 0.1 -> 10, 0.01 -> 20 etc.)
    phred_t toPhredQ( const double probErr )
    {
        return std::min( phredCoefficient * log10( probErr ), constants::maxPhredScore );
    }

    unsigned int roundPhred( const double phred ) { return static_cast< unsigned int >( std::lround( phred ) ); }

    // Converts a "PHRED" quality score to an error probability
    double fromPhredQ( const double phredQ ) { return pow( 10.0, phredQ / phredCoefficient ); }

    double variantSupportPerRead( const double prior, const double posterior, const int64_t variantSupportCount )
    {
        if ( variantSupportCount == 0L or ( posterior + prior ) >= constants::maxPhredScore )
        {
            return std::numeric_limits< double >::quiet_NaN();
        }
        else
        {
            return ( posterior + prior ) / static_cast< double >( variantSupportCount );
        }
    }
}
}
