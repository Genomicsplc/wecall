// All content Copyright (C) 2018 Genomics plc
#ifndef COMMON_TYPES_HPP
#define COMMON_TYPES_HPP

#include <iostream>
#include <iterator>
#include <sstream>

namespace echidna
{
namespace corrector
{
    const static double phred_to_p_cache[100] = {1.0,
                                                 0.794328234724,
                                                 0.63095734448,
                                                 0.501187233627,
                                                 0.398107170553,
                                                 0.316227766017,
                                                 0.251188643151,
                                                 0.199526231497,
                                                 0.158489319246,
                                                 0.125892541179,
                                                 0.1,
                                                 0.0794328234724,
                                                 0.063095734448,
                                                 0.0501187233627,
                                                 0.0398107170553,
                                                 0.0316227766017,
                                                 0.0251188643151,
                                                 0.0199526231497,
                                                 0.0158489319246,
                                                 0.0125892541179,
                                                 0.01,
                                                 0.00794328234724,
                                                 0.0063095734448,
                                                 0.00501187233627,
                                                 0.00398107170553,
                                                 0.00316227766017,
                                                 0.00251188643151,
                                                 0.00199526231497,
                                                 0.00158489319246,
                                                 0.00125892541179,
                                                 0.001,
                                                 0.000794328234724,
                                                 0.00063095734448,
                                                 0.000501187233627,
                                                 0.000398107170553,
                                                 0.000316227766017,
                                                 0.000251188643151,
                                                 0.000199526231497,
                                                 0.000158489319246,
                                                 0.000125892541179,
                                                 0.0001,
                                                 7.94328234724e-05,
                                                 6.3095734448e-05,
                                                 5.01187233627e-05,
                                                 3.98107170553e-05,
                                                 3.16227766017e-05,
                                                 2.51188643151e-05,
                                                 1.99526231497e-05,
                                                 1.58489319246e-05,
                                                 1.25892541179e-05,
                                                 1e-05,
                                                 7.94328234724e-06,
                                                 6.3095734448e-06,
                                                 5.01187233627e-06,
                                                 3.98107170553e-06,
                                                 3.16227766017e-06,
                                                 2.51188643151e-06,
                                                 1.99526231497e-06,
                                                 1.58489319246e-06,
                                                 1.25892541179e-06,
                                                 1e-06,
                                                 7.94328234724e-07,
                                                 6.3095734448e-07,
                                                 5.01187233627e-07,
                                                 3.98107170553e-07,
                                                 3.16227766017e-07,
                                                 2.51188643151e-07,
                                                 1.99526231497e-07,
                                                 1.58489319246e-07,
                                                 1.25892541179e-07,
                                                 1e-07,
                                                 7.94328234724e-08,
                                                 6.3095734448e-08,
                                                 5.01187233627e-08,
                                                 3.98107170553e-08,
                                                 3.16227766017e-08,
                                                 2.51188643151e-08,
                                                 1.99526231497e-08,
                                                 1.58489319246e-08,
                                                 1.25892541179e-08,
                                                 1e-08,
                                                 7.94328234724e-09,
                                                 6.3095734448e-09,
                                                 5.01187233627e-09,
                                                 3.98107170553e-09,
                                                 3.16227766017e-09,
                                                 2.51188643151e-09,
                                                 1.99526231497e-09,
                                                 1.58489319246e-09,
                                                 1.25892541179e-09,
                                                 1e-09,
                                                 7.94328234724e-10,
                                                 6.3095734448e-10,
                                                 5.01187233627e-10,
                                                 3.98107170553e-10,
                                                 3.16227766017e-10,
                                                 2.51188643151e-10,
                                                 1.99526231497e-10,
                                                 1.58489319246e-10,
                                                 1.25892541179e-10};

    const std::size_t kmerSize = 7;
    const std::size_t padding = 1;
    const double MINIMUM_KMER_PRIOR = 2e-3;  ///< kmers with setPrior less than this are not considered

    template < int k >
    using kmer_t = std::array< char, k >;

    template < int k >
    struct kmerhash_t
    {
        std::size_t operator()( const kmer_t< k > & kmer ) const
        {
            uint32_t hash = 0;
            for ( int i = 0; i < k; i++ )
            {
                // inspired by stackoverflow.com/questions/3062746
                hash += 12345 + ( ( kmer[i] & 6 ) >> 1 );
                hash *= 1103515245;
            }
            return hash;
        }
    };

    template < int k, int padding >
    using extKmer_t = std::array< char, k + 2 * padding >;

    struct ErrorCountData
    {
        double errorOpportunity;
        double errorCount;
    };

    // http://stackoverflow.com/questions/16782746/what-is-faster-than-stdpow
    inline double fastPow( const double a, const double b )
    {
        union
        {
            double d;
            int x[2];
        } u = {a};

        u.x[1] = (int)( b * ( u.x[1] - 1072632447 ) + 1072632447 );
        u.x[0] = 0;

        return u.d;
    }

    inline double phred_to_p( const int phred )
    {
        return phred_to_p_cache[phred];
        // return std::pow( 0.1, phred / 10.0 );
    }

    template < typename S, typename T >
    void show2( S s, T a )
    {
        std::cout << s;
        std::copy( a.begin(), a.end(), std::ostream_iterator< typename T::value_type >( std::cout, " " ) );
        std::cout << std::endl;
    }

    template < typename T >
    void show( T a )
    {
        std::copy( a.begin(), a.end(), std::ostream_iterator< typename T::value_type >( std::cout, "" ) );
    }

    template < typename T >
    std::string show_string( T a )
    {
        std::stringstream result;
        for ( auto c_it = a.cbegin(); c_it != a.cend(); ++c_it )
        {
            result << *c_it;
        }
        return result.str();
    }

}  // namespace corrector
}  // namespace echidna

#endif
