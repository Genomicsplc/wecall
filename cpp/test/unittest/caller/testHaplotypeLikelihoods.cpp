#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "utils/matrix.hpp"
#include "caller/haplotypeLikelihoods.hpp"

using echidna::utils::matrix_t;

echidna::utils::matrix_t getHaplotypeMatrixFromVecOfVecs( const std::vector< std::vector< double > > & values )
{
    assert( values.size() > 0 );

    echidna::utils::matrix_t matrix_t( values.size(), values[0].size() );

    for ( std::size_t rowIndex = 0; rowIndex < matrix_t.size1(); ++rowIndex )
    {
        const auto & row = values[rowIndex];
        assert( row.size() == matrix_t.size2() );

        std::copy( row.begin(), row.end(), matrix_t.data().begin() + matrix_t.size2() * rowIndex );
    }
    return matrix_t;
}

BOOST_AUTO_TEST_CASE( shouldComputeHaplotypeFrequencies )
{
    std::vector< std::vector< double > > values = {{0.5, 0.3, 0.2}, {0.2, 0.7, 0.1}};
    matrix_t haplotypeLikelihoods = getHaplotypeMatrixFromVecOfVecs( values );
    assert( haplotypeLikelihoods.size1() == 2 );
    assert( haplotypeLikelihoods.size2() == 3 );

    auto result = echidna::caller::computeHaplotypeFrequencies( haplotypeLikelihoods, {} );

    BOOST_CHECK_EQUAL( result.size(), haplotypeLikelihoods.size2() );
    BOOST_CHECK_CLOSE( result[0], 0.7, 1e-6 );
    BOOST_CHECK_CLOSE( result[1], 1.0, 1e-6 );
    BOOST_CHECK_CLOSE( result[2], 0.3, 1e-6 );
}

BOOST_AUTO_TEST_CASE( shouldComputeHaplotypeFrequenciesWhenNoDataConsidered )
{
    std::vector< std::vector< double > > values = {{0.5, 0.3}, {0.2, 0.7}};
    matrix_t haplotypeLikelihoods = getHaplotypeMatrixFromVecOfVecs( values );
    assert( haplotypeLikelihoods.size1() == 2 );
    assert( haplotypeLikelihoods.size2() == 2 );

    auto result = echidna::caller::computeHaplotypeFrequencies( haplotypeLikelihoods, {0, 1} );

    BOOST_CHECK_EQUAL( result.size(), 2 );
    BOOST_CHECK_CLOSE( result[0], 0.0, 1e-6 );
    BOOST_CHECK_CLOSE( result[1], 0.0, 1e-6 );
}

BOOST_AUTO_TEST_CASE( shouldComputeHaplotypeFrequenciesExcludingSpecifiedHaplotypes )
{
    std::vector< std::vector< double > > values = {{0.5, 0.3, 0.2}, {0.2, 0.7, 0.1}};
    matrix_t haplotypeLikelihoods = getHaplotypeMatrixFromVecOfVecs( values );
    assert( haplotypeLikelihoods.size1() == 2 );
    assert( haplotypeLikelihoods.size2() == 3 );

    auto result = echidna::caller::computeHaplotypeFrequencies( haplotypeLikelihoods, {1} );

    BOOST_CHECK_EQUAL( result.size(), haplotypeLikelihoods.size2() );
    BOOST_CHECK_CLOSE( result[0], 0.5 / 0.7 + 0.2 / 0.3, 1e-6 );
    BOOST_CHECK_CLOSE( result[1], 0.0, 1e-6 );
    BOOST_CHECK_CLOSE( result[2], 0.2 / 0.7 + 0.1 / 0.3, 1e-6 );
}
