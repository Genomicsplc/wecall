// All content Copyright (C) 2018 Genomics plc
#include "variant/variantFilter.hpp"
#include "variant/variantContainer.hpp"
#include "common.hpp"
#include "utils/logging.hpp"

namespace wecall
{
namespace variant
{
    //-----------------------------------------------------------------------------------------

    variantSet_t VariantFilter::getSortedFilteredVariants( const caller::Region & blockRegion,
                                                           const VariantContainer & varContainer ) const
    {
        const auto variants = varContainer.getVariants();

        variantSet_t filteredVars;
        for ( const auto varPtr : variants )
        {
            if ( blockRegion.contains( varPtr->region() ) and variantPassesFilters( varPtr, varContainer ) )
            {
                filteredVars.insert( varPtr );
            }
        }

        return filteredVars;
    }

    //-----------------------------------------------------------------------------------------

    bool VariantFilter::variantHasMinReadsAcrossSamples( varPtr_t varPtr, const VariantContainer & varContainer ) const
    {
        const auto nSupportingReads = varContainer.totalReadsSupportingVariant( varPtr );
        if ( false )
        {
            ECHIDNA_LOG( SUPER_DEBUG, "var " << varPtr->toString() << " has " << nSupportingReads
                                             << " supporting reads" );
        }
        return m_minReads <= nSupportingReads;
    }

    //-----------------------------------------------------------------------------------------

    bool VariantFilter::variantHasSampleWithMinPercentReadCoverage( varPtr_t varPtr,
                                                                    const VariantContainer & varContainer ) const
    {
        const auto percentSupportingReads = varContainer.maxReadPercentVariantCoverage( varPtr );
        if ( false )
        {
            ECHIDNA_LOG( SUPER_DEBUG, "var " << varPtr->toString() << " has max % support of "
                                             << percentSupportingReads );
        }
        return m_minPerSamplePercentage <= percentSupportingReads;
    }

    //-----------------------------------------------------------------------------------------

    bool VariantFilter::variantPassesFilters( varPtr_t varPtr, const VariantContainer & varContainer ) const
    {
        if ( varPtr->neverFilter() )
        {
            return true;
        }

        if ( not variantHasMinReadsAcrossSamples( varPtr, varContainer ) )
        {
            return false;
        }

        return variantHasSampleWithMinPercentReadCoverage( varPtr, varContainer );
    }

    //-----------------------------------------------------------------------------------------
}
}
