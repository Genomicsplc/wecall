// All content Copyright (C) 2018 Genomics plc
#include "utils/multinomialCoefficients.hpp"

#include <numeric>

unsigned int factorial( unsigned int n )
{
    unsigned int result = 1;
    while ( n )
    {
        result *= n--;
    }
    return result;
}

unsigned int multinomial_coefficient( std::vector< unsigned int > k )
{
    const unsigned int n = std::accumulate( k.cbegin(), k.cend(), 0u );

    unsigned int den = 1;
    for ( const auto & k_i : k )
    {
        den *= factorial( k_i );
    }

    return factorial( n ) / den;
}
