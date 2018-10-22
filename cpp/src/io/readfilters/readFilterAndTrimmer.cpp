// All content Copyright (C) 2018 Genomics plc
#include "readFilterAndTrimmer.hpp"

#include "common.hpp"
#include "utils/logging.hpp"

#include "io/read.hpp"
#include "io/readfilters/rangeFilter.hpp"
#include "io/readfilters/booleanFilter.hpp"
#include "io/readfilters/mapQualityFilter.hpp"
#include "io/readfilters/baseQualityFilter.hpp"
#include "io/readfilters/shortReadFilter.hpp"

namespace wecall
{
namespace io
{
    ReadFilterAndTrimmer::ReadFilterAndTrimmer( const caller::params::Filters & filterParams )
        : m_overlapTrim( filterParams.m_overlapTrim ),
          m_shortReadTrim( filterParams.m_shortReadTrim ),
          m_noSimilarReads( filterParams.m_noSimilarReadsFilter )
    {
        ECHIDNA_LOG( DEBUG, "Initialising read level filters in ReadFiltersManager" );

        m_filters.emplace_back( std::make_shared< BooleanFilter >( &Read::isUnMapped, "isUnmapped", true ) );

        m_filters.emplace_back( std::make_shared< BooleanFilter >( &Read::isSecondary, "isSecondary", true ) );

        // Configurable quality filters
        m_filters.emplace_back(
            std::make_shared< BaseQualFilter >( filterParams.m_baseCallFilterQ, filterParams.m_baseCallFilterN ) );

        // Configurable boolean filters
        if ( filterParams.m_duplicatesFilter )
        {
            m_filters.emplace_back( std::make_shared< BooleanFilter >( &Read::isDuplicate, "isDuplicate", true ) );
        }

        if ( not filterParams.m_allowImproperPairs )
        {
            m_filters.emplace_back( std::make_shared< BooleanFilter >( &Read::isProperPair, "isProperPair" ) );
        }

        if ( filterParams.m_noMatesFilter )
        {
            m_filters.emplace_back(
                std::make_shared< BooleanFilter >( &Read::isMateUnMapped, "isMateUnmapped", true ) );
        }

        if ( filterParams.m_shortReadFilter )
        {
            m_filters.emplace_back( std::make_shared< ShortReadFilter >() );
        }
    }

    bool ReadFilterAndTrimmer::isSimilarToPrevious( readPtr_t read ) const
    {
        bool isSimilarToPrevious = false;
        if ( m_previousRead != nullptr )
        {
            // AR: Two read-pairs are duplicates if their fragments span the same locus.
            if ( m_previousRead->getStartPos() == read->getStartPos() and
                 m_previousRead->getInsertSize() == read->getInsertSize() )
            {
                isSimilarToPrevious = true;
            }
        }
        m_previousRead = read;
        return isSimilarToPrevious;
    }

    bool ReadFilterAndTrimmer::trimAndFilter( readPtr_t read ) const
    {
        trim( read );
        if ( passesFilters( read ) and hasLength( read ) )
        {
            if ( m_noSimilarReads )
            {
                return not this->isSimilarToPrevious( read );
            }
            else
            {
                return true;
            }
        }
        else
        {
            return false;
        }
    }

    bool ReadFilterAndTrimmer::hasLength( readPtr_t read ) const
    {
        return ( read->getStartPos() < read->getAlignedEndPos() );
    }

    void ReadFilterAndTrimmer::trim( readPtr_t read ) const
    {
        if ( m_overlapTrim )
        {
            read->trimOverlap();
        }

        if ( m_shortReadTrim )
        {
            read->trimReadOfShortFragment();
        }
    }

    bool ReadFilterAndTrimmer::passesFilters( readPtr_t read ) const
    {
        bool passesAll = true;
        for ( auto & filter : m_filters )
        {
            // Note:    Test all filters even if one has already failed, so that counts are
            //            complete and not affected by filter test order.

            // TOD0 - change passes filter to take either a pointer or a ref to a pointer
            if ( not filter->passesFilter( *read ) )
            {
                passesAll = false;
            }
        }

        return passesAll;
    }
}
}
