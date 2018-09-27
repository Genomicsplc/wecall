// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <string>
#include <bits/shared_ptr.h>

#include "variant/genotype.hpp"
#include "variant/haplotype.hpp"
#include "variant/type/variant.hpp"

using namespace echidna::variant;
using echidna::variant::variantSet_t;
using echidna::utils::ReferenceSequence;
using echidna::caller::Region;
using echidna::variant::HaplotypeVector;

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( haploidGenotype )
{
    using namespace echidna::variant;

    const std::string chrom( "1" );

    const auto reference = std::make_shared< ReferenceSequence >( Region( chrom, 1000010, 1000011 ), "T" );

    const varPtr_t var1 = std::make_shared< Variant >( reference, reference->region(), "A" );
    const varPtr_t var2 = std::make_shared< Variant >( reference, reference->region(), "G" );

    const variantSet_t variants1 = {var1};

    HaplotypeVector haplotypes( reference->region(), reference );
    haplotypes.push_back( variants1 );

    const Genotype theGenotype( haplotypeAndCount_t{{0, 1}} );

    BOOST_CHECK_EQUAL( theGenotype.nStrandsContainingVariant( haplotypes, var1 ), 1 );
    BOOST_CHECK_EQUAL( theGenotype.nStrandsContainingVariant( haplotypes, var2 ), 0 );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( diploidGenotype )
{
    using namespace echidna::variant;

    const std::string chrom( "1" );

    const auto reference = std::make_shared< ReferenceSequence >( Region( chrom, 1000010, 1000012 ), "TT" );
    const varPtr_t var1 = std::make_shared< Variant >( reference, Region( chrom, 1000010, 1000011 ), "A" );
    const varPtr_t var2 = std::make_shared< Variant >( reference, Region( chrom, 1000011, 1000012 ), "G" );

    const variantSet_t variants1 = {var1};
    const variantSet_t variants2 = {var1, var2};

    HaplotypeVector haplotypes( reference->region(), reference );
    haplotypes.push_back( variants1 );
    haplotypes.push_back( variants2 );

    const Genotype theGenotype( haplotypeAndCount_t{{0, 1}, {1, 1}} );

    BOOST_CHECK_EQUAL( theGenotype.nStrandsContainingVariant( haplotypes, var1 ), 2 );
    BOOST_CHECK_EQUAL( theGenotype.nStrandsContainingVariant( haplotypes, var2 ), 1 );
}

//-------------------------------------------------------------------------------------------------
