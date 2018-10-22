// All content Copyright (C) 2018 Genomics plc
#include <cassert>

#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include "caller/region.hpp"
#include "io/fastaFile.hpp"
#include "utils/exceptions.hpp"
#include "utils/logging.hpp"
#include "utils/interval.hpp"
#include "common.hpp"

namespace wecall
{
namespace caller
{

    void Region::combine( const Region & other )
    {
        WECALL_ASSERT( ( m_contig == other.m_contig ),
                        "Can not combine regions " + this->toString() + " & " + other.toString() );
        m_interval.combine( other.m_interval );
    }

    bool Region::contains( const Region & other ) const
    {
        return m_contig == other.m_contig and m_interval.contains( other.m_interval );
    }

    bool Region::contains( const std::string & otherContig, int64_t otherPos ) const
    {
        return m_contig == otherContig and m_interval.contains( otherPos );
    }

    bool Region::overlapsOrTouches( const Region & other ) const
    {
        if ( m_contig != other.m_contig )
        {
            return false;
        }

        return m_interval.overlapsOrTouches( other.m_interval );
    }

    Region Region::getIntersect( const Region & other ) const
    {
        WECALL_ASSERT( this->contig() == other.contig(), "Can't intersect regions from different contigs" );

        const auto newInterval = this->interval().getIntersect( other.interval() );
        return Region( this->contig(), newInterval.start(), newInterval.end() );
    }

    std::string Region::toString() const
    {
        std::stringstream stream;
        stream << *this;
        return stream.str();
    }

    bool Region::operator<( const Region & other ) const
    {
        if ( m_contig != other.m_contig )
        {
            return m_contig < other.m_contig;
        }
        return m_interval < other.m_interval;
    }

    bool Region::operator==( const Region & other ) const
    {
        return m_contig == other.m_contig and m_interval == other.m_interval;
    }

    bool Region::operator!=( const Region & other ) const { return not( *this == other ); }

    Region SetRegions::getSpan() const
    {
        return Region( m_regions.cbegin()->contig(), m_regions.cbegin()->start(), m_regions.rbegin()->end() );
    }

    std::string SetRegions::toString() const
    {
        std::stringstream sstr;
        sstr << "Regions: {";
        for ( const auto & region : m_regions )
        {
            sstr << region << ",";
        }
        sstr << "}";
        return sstr.str();
    }

    bool SetRegions::operator==( const SetRegions & other ) const
    {
        return other.size() == this->size() and
               std::equal( m_regions.begin(), m_regions.end(), other.m_regions.begin() );
    }

    bool SetRegions::operator<( const SetRegions & other ) const
    {
        if ( m_regions.size() != other.m_regions.size() )
        {
            return m_regions.size() < other.m_regions.size();
        }
        for ( auto it1 = m_regions.begin(), it2 = other.m_regions.begin();
              it1 != m_regions.end() and it2 != other.m_regions.end(); ++it1, ++it2 )
        {
            if ( *it1 != *it2 )
            {
                return *it1 < *it2;
            }
        }
        return false;
    }

    bool SetRegions::allSameContig() const
    {
        if ( m_regions.size() <= 1 )
        {
            return true;
        }
        const auto firstContig = m_regions.cbegin()->contig();
        const auto sameContig = [firstContig]( const Region & region ) -> bool
        {
            return region.contig() == firstContig;
        };
        return std::all_of( m_regions.cbegin(), m_regions.cend(), sameContig );
    }

    void SetRegions::fill( const int64_t fillDistance )
    {
        if ( m_regions.size() < 2 )
        {
            return;
        }
        std::vector< Region > newRegions;
        for ( auto it1 = m_regions.begin(), it2 = std::next( it1 ); it2 != m_regions.end(); ++it1, ++it2 )
        {
            if ( it1->contig() == it2->contig() and it2->start() - it1->end() <= fillDistance )
            {
                newRegions.emplace_back( it1->contig(), it1->start(), it2->end() );
            }
        }
        for ( const auto & reg : newRegions )
        {
            this->insert( reg );
        }
    }

    void SetRegions::insert( const Region & paddedRegion )
    {
        const auto regionOverlaps = [paddedRegion]( const caller::Region & other )
        {
            return paddedRegion.overlapsOrTouches( other );
        };
        const auto contig = paddedRegion.contig();
        auto minPos = paddedRegion.start();
        auto maxPos = paddedRegion.end();
        for ( auto it = m_regions.begin(); it != m_regions.end(); )
        {
            if ( regionOverlaps( *it ) )
            {
                minPos = std::min( it->start(), minPos );
                maxPos = std::max( it->end(), maxPos );
                m_regions.erase( it++ );
            }
            else
            {
                ++it;
            }
        }
        m_regions.insert( caller::Region( contig, minPos, maxPos ) );
    }

    std::ostream & operator<<( std::ostream & out, const Region & theRegion )
    {
        if ( theRegion.hasNoRange() )
        {
            return out << theRegion.contig();
        }
        // else
        // This format is required for BAM via Samtools.
        return out << theRegion.contig() << ":" << theRegion.start() << "-" << theRegion.end();
    }

    std::ostream & operator<<( std::ostream & out, const SetRegions & theRegions )
    {
        // This format is required for BAM via Samtools.
        return out << theRegions.toString();
    }
}
}
