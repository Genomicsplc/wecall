// All content Copyright (C) 2018 Genomics plc
#include "caller/region.hpp"
#include "variant/variantContainer.hpp"
#include "variant/variantFilter.hpp"
#include "utils/referenceSequence.hpp"
#include "utils/sequence.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

using echidna::variant::Variant;
using echidna::variant::VariantFilter;
using echidna::variant::VariantContainer;
using echidna::utils::ReferenceSequence;
using echidna::utils::BasePairSequence;
using echidna::caller::Region;

std::shared_ptr< Variant > unfilterable_variant( const std::string & chrom,
                                                 const int64_t start,
                                                 const int64_t end,
                                                 const std::string & ref,
                                                 const std::string & alt )
{
    auto referenceSequence =
        std::make_shared< ReferenceSequence >( Region( chrom, start, end ), BasePairSequence( ref ) );
    std::shared_ptr< Variant > var =
        std::make_shared< Variant >( referenceSequence, referenceSequence->region(), BasePairSequence( alt ) );
    var->neverFilter();
    return var;
}

BOOST_AUTO_TEST_CASE( testVariantFilterAllowsOnlyVariantsInGivenRegion )
{
    VariantContainer container( 0, 0 );

    container.addCandidateVariant( unfilterable_variant( "1", 95, 96, "C", "T" ), 0.0 );
    container.addCandidateVariant( unfilterable_variant( "1", 100, 101, "C", "T" ), 0.0 );
    container.addCandidateVariant( unfilterable_variant( "1", 105, 106, "C", "T" ), 0.0 );
    container.addCandidateVariant( unfilterable_variant( "1", 110, 111, "C", "T" ), 0.0 );

    const auto filter = VariantFilter( 2, 2 );
    const auto actual = filter.getSortedFilteredVariants( Region( "1", 100, 106 ), container );

    BOOST_CHECK_EQUAL( 2, actual.size() );
}
