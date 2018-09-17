#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <string>
#include <bits/shared_ptr.h>

#include "variant/genotype.hpp"
#include "variant/haplotype.hpp"
#include "variant/type/variant.hpp"
#include "caller/diploid/diploidAnnotate.hpp"

using echidna::variant::GenotypeVector;
using echidna::variant::HaplotypeVector;
using echidna::variant::genotypePtr_t;
using echidna::variant::varPtr_t;
using echidna::variant::variantSet_t;
using echidna::variant::Haplotype;
using echidna::utils::ReferenceSequence;
using echidna::caller::Region;
using echidna::variant::Variant;

BOOST_AUTO_TEST_CASE( testConstructsAllPloidy3CombinationsFromHaplotypeVector )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );

    varPtr_t snp1 = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );
    varPtr_t snp2 = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "G" );
    varPtr_t snp3 = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "T" );
    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );
    haplotypeVector.push_back( variantSet_t{} );
    haplotypeVector.push_back( variantSet_t{snp1} );
    haplotypeVector.push_back( variantSet_t{snp2} );
    haplotypeVector.push_back( variantSet_t{snp3} );

    const auto ploidy = 3;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1, 2, 3}, {} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 20 );  // multinomialElements(4, 3)
    BOOST_CHECK_EQUAL( genotypeVector.ploidy(), ploidy );

    // Check genotype uniqueness
    std::set< std::vector< std::size_t > > haplotypeIndicesSet;
    for ( std::size_t i = 0; i < genotypeVector.size(); ++i )
    {
        const auto haplotypeIndices = genotypeVector.getHaplotypeIndices( i );
        BOOST_REQUIRE_EQUAL( haplotypeIndices.size(), ploidy );

        haplotypeIndicesSet.insert( haplotypeIndices );
    }
    BOOST_CHECK_EQUAL( haplotypeIndicesSet.size(), 20 );
}

BOOST_AUTO_TEST_CASE( testPloidy3GenotypeCombinations )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );

    varPtr_t snp = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );
    HaplotypeVector haplotypeVector( Region( "1", 0, 1 ), referenceSequence );
    haplotypeVector.push_back( variantSet_t{} );
    haplotypeVector.push_back( variantSet_t{snp} );

    const auto ploidy = 3;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1}, {} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 4 );

    BOOST_CHECK_EQUAL( genotypeVector[0]->nCombinationsThisGenotype(), 1 );
    BOOST_CHECK_EQUAL( genotypeVector[1]->nCombinationsThisGenotype(), 3 );
    BOOST_CHECK_EQUAL( genotypeVector[2]->nCombinationsThisGenotype(), 3 );
    BOOST_CHECK_EQUAL( genotypeVector[3]->nCombinationsThisGenotype(), 1 );
}

BOOST_AUTO_TEST_CASE( testPloidy3GenotypeVectorOrdering )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 2 ), "AA" );

    const auto snp0 = std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "C" );

    HaplotypeVector haplotypeVector( Region( "1", 0, 2 ), referenceSequence );
    auto regionReference = haplotypeVector.regionReferenceSequence();
    haplotypeVector.push_back( variantSet_t{} );      //  hap0
    haplotypeVector.push_back( variantSet_t{snp0} );  //  hap1

    const auto ploidy = 3;

    const GenotypeVector genotypeVector( ploidy, haplotypeVector, {0, 1}, {} );

    BOOST_REQUIRE_EQUAL( genotypeVector.size(), 4 );

    // Ref-ref-ref
    BOOST_REQUIRE_EQUAL( genotypeVector[0]->nStrandsContainingVariant( haplotypeVector, snp0 ), 0 );

    // Ref-ref-hap1
    BOOST_REQUIRE_EQUAL( genotypeVector[1]->nStrandsContainingVariant( haplotypeVector, snp0 ), 1 );

    // Ref-hap1-hap1
    BOOST_REQUIRE_EQUAL( genotypeVector[2]->nStrandsContainingVariant( haplotypeVector, snp0 ), 2 );

    // hap1-hap1-hap1
    BOOST_REQUIRE_EQUAL( genotypeVector[3]->nStrandsContainingVariant( haplotypeVector, snp0 ), 3 );
}
