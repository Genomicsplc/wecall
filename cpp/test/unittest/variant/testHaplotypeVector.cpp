#include "variant/haplotype.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <string>

#include "variant/type/variant.hpp"
#include "utils/referenceSequence.hpp"
#include "common.hpp"

using echidna::utils::ReferenceSequence;
using echidna::utils::referenceSequencePtr_t;
using echidna::caller::Region;
using echidna::variant::HaplotypeVector;
using echidna::variant::varPtr_t;
using echidna::variant::Haplotype;
using echidna::variant::Variant;
using echidna::variant::variantSet_t;

BOOST_AUTO_TEST_CASE( constructorShouldRaiseExceptionIfRegionsNotCompatible )
{

    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );
    const Region incompatibleRegion( "2", 0, 1 );
    BOOST_CHECK_THROW( HaplotypeVector( incompatibleRegion, referenceSequence ), echidna::utils::echidna_exception );
}

BOOST_AUTO_TEST_CASE( testShouldBeAbleToAddHaplotypeForSNP )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );
    HaplotypeVector haplotypeVector( referenceSequence->region(), referenceSequence );

    BOOST_CHECK( haplotypeVector.empty() );

    haplotypeVector.push_back( {std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "T" )} );

    BOOST_CHECK( not haplotypeVector.empty() );
    BOOST_CHECK_EQUAL( haplotypeVector.size(), 1 );
    BOOST_CHECK_EQUAL( haplotypeVector[0].paddedSequences().front(), "T" );
}

BOOST_AUTO_TEST_CASE( testShouldRaiseForOverlappingVariants )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );
    HaplotypeVector haplotypeVector( referenceSequence->region(), referenceSequence );

    BOOST_CHECK( haplotypeVector.empty() );

    BOOST_CHECK_THROW( haplotypeVector.push_back(
                           {std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "T" ),
                            std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "G" )} ),
                       echidna::utils::echidna_exception );
}

BOOST_AUTO_TEST_CASE( testShouldSortByAltString )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );
    HaplotypeVector haplotypeVector( referenceSequence->region(), referenceSequence );
    haplotypeVector.push_back( {std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "T" )} );
    haplotypeVector.push_back( {std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "G" )} );

    BOOST_CHECK_EQUAL( haplotypeVector[0].paddedSequences().front(), "T" );
    BOOST_CHECK_EQUAL( haplotypeVector[1].paddedSequences().front(), "G" );

    haplotypeVector.sort();

    BOOST_CHECK_EQUAL( haplotypeVector[0].paddedSequences().front(), "G" );
    BOOST_CHECK_EQUAL( haplotypeVector[1].paddedSequences().front(), "T" );
}

BOOST_AUTO_TEST_CASE( testShouldPadHaplotypeWithReferenceString )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 2 ), "AA" );
    HaplotypeVector haplotypeVector( referenceSequence->region(), referenceSequence );

    haplotypeVector.push_back( {std::make_shared< Variant >( referenceSequence, Region( "1", 1, 1 ), "T" )} );

    BOOST_CHECK_EQUAL( haplotypeVector[0].paddedSequences().front(), "ATA" );
}

BOOST_AUTO_TEST_CASE( testShouldFailForOutOfRangePaddedSequenceIndex )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 2 ), "AA" );
    HaplotypeVector haplotypeVector( referenceSequence->region(), referenceSequence );

    BOOST_CHECK_THROW( haplotypeVector[0].paddedSequences().front(), std::out_of_range );
}

BOOST_AUTO_TEST_CASE( testShouldFailForOutOfRangeHaplotypeIndex )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 2 ), "AA" );
    HaplotypeVector haplotypeVector( referenceSequence->region(), referenceSequence );

    BOOST_CHECK_THROW( haplotypeVector[0], std::out_of_range );
}

BOOST_AUTO_TEST_CASE( testHaplotypeReferencePadding )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "0123456789" );

    HaplotypeVector haplotypeVector( Region( "1", 4, 6 ), referenceSequence );

    haplotypeVector.push_back( variantSet_t{} );

    BOOST_CHECK_EQUAL( referenceSequence->sequence(), haplotypeVector[0].paddedSequences().front() );
}

BOOST_AUTO_TEST_CASE( testHaplotypeReferencePaddingWithSNP )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 10 ), "AAAAAAAAAA" );

    HaplotypeVector haplotypeVector( Region( "1", 4, 6 ), referenceSequence );

    haplotypeVector.push_back( {std::make_shared< Variant >( referenceSequence, Region( "1", 4, 5 ), "G" )} );

    BOOST_CHECK_EQUAL( "AAAAGAAAAA", haplotypeVector[0].paddedSequences().front() );
}

BOOST_AUTO_TEST_CASE( testShouldBeAbleToIterateOverAndIndexIntoVector )
{
    auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 1 ), "A" );
    HaplotypeVector haplotypeVector( referenceSequence->region(), referenceSequence );
    haplotypeVector.push_back( {std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "G" )} );
    haplotypeVector.push_back( {std::make_shared< Variant >( referenceSequence, referenceSequence->region(), "T" )} );

    auto it = haplotypeVector.begin();
    BOOST_CHECK( it != haplotypeVector.end() );
    BOOST_CHECK_EQUAL( *it, haplotypeVector[0] );

    ++it;
    BOOST_CHECK( it != haplotypeVector.end() );
    BOOST_CHECK_EQUAL( *it, haplotypeVector[1] );

    ++it;
    BOOST_CHECK( it == haplotypeVector.end() );
}

BOOST_AUTO_TEST_CASE( testShouldMergeMatchingInsertionsTogether )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 2 ), "AA" );
    HaplotypeVector haplotypeVector( referenceSequence->region(), referenceSequence );

    haplotypeVector.push_back( {std::make_shared< Variant >( referenceSequence, Region( "1", 0, 0 ), "A" )} );
    haplotypeVector.push_back( {std::make_shared< Variant >( referenceSequence, Region( "1", 1, 1 ), "A" )} );

    haplotypeVector.sort();
    haplotypeVector.merge();

    BOOST_CHECK_EQUAL( haplotypeVector.size(), 1 );

    BOOST_CHECK( haplotypeVector[0].containsVariant(
        {std::make_shared< Variant >( referenceSequence, Region( "1", 0, 0 ), "A" )} ) );
    BOOST_CHECK( not haplotypeVector[0].containsVariant(
        {std::make_shared< Variant >( referenceSequence, Region( "1", 1, 1 ), "A" )} ) );
}

BOOST_AUTO_TEST_CASE( testShouldMergeMatchingDeletionsTogether )
{
    const auto referenceSequence = std::make_shared< ReferenceSequence >( Region( "1", 0, 2 ), "AA" );
    HaplotypeVector haplotypeVector( referenceSequence->region(), referenceSequence );

    haplotypeVector.push_back( {std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "" )} );
    haplotypeVector.push_back( {std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "" )} );

    haplotypeVector.sort();
    haplotypeVector.merge();

    BOOST_CHECK_EQUAL( haplotypeVector.size(), 1 );

    BOOST_CHECK( haplotypeVector[0].containsVariant(
        {std::make_shared< Variant >( referenceSequence, Region( "1", 0, 1 ), "" )} ) );
    BOOST_CHECK( not haplotypeVector[0].containsVariant(
        {std::make_shared< Variant >( referenceSequence, Region( "1", 1, 2 ), "" )} ) );
}
