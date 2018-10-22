// All content Copyright (C) 2018 Genomics plc
#ifndef CALLER_REGION_HPP
#define CALLER_REGION_HPP

#include <string>
#include <cstdint>
#include <sstream>
#include <vector>
#include <iostream>
#include <map>
#include <set>

#include "utils/interval.hpp"

namespace wecall
{

namespace caller
{
    const std::string chromMatch = "([\\w,\\|,\\.]+)";
    const std::string posMatch = "(\\d+)";

    class Region
    {
    public:
        Region( std::string contigIn, int64_t startIn, int64_t endIn )
            : m_contig( contigIn ), m_interval( startIn, endIn )
        {
        }

        Region( std::string contigIn, utils::Interval interval ) : m_contig( contigIn ), m_interval( interval ) {}

        std::string contig() const { return m_contig; }
        utils::Interval interval() const { return m_interval; }
        int64_t start() const { return m_interval.start(); }
        int64_t end() const { return m_interval.end(); }

        bool hasNoRange() const { return ( end() == 0 ); }
        int64_t size() const { return m_interval.size(); }

        bool contains( const Region & other ) const;
        bool contains( const std::string & otherContig, int64_t otherPos ) const;
        bool overlaps( const Region & other ) const
        {
            return m_contig == other.m_contig and m_interval.overlaps( other.m_interval );
        }
        bool overlapsOrTouches( const Region & other ) const;

        // Take maximal span to two regions
        void combine( const Region & other );

        bool operator<( const Region & other ) const;
        bool operator==( const Region & other ) const;
        bool operator!=( const Region & other ) const;

        Region getPadded( const int64_t & paddingSize ) const
        {
            return Region( m_contig, m_interval.getPadded( paddingSize ) );
        }

        Region getIntersect( const Region & other ) const;

        // Utility function so we can std::cout a region
        friend std::ostream & operator<<( std::ostream & out, const Region & theRegion );

        std::string toString() const;

    private:
        std::string m_contig;
        utils::Interval m_interval;
    };

    class SetRegions
    {
    public:
        SetRegions() {}
        SetRegions( const Region & region ) : m_regions( {region} ) {}

        bool operator==( const SetRegions & other ) const;
        bool operator!=( const SetRegions & other ) const { return not( *this == other ); }
        bool operator<( const SetRegions & other ) const;

        bool allSameContig() const;
        void insert( const Region & region );
        void fill( const int64_t fillDistance );
        Region getSpan() const;

        bool overlaps( const caller::Region & region ) const
        {
            for ( auto & r : m_regions )
            {
                if ( r.overlaps( region ) )
                {
                    return true;
                }
            }
            return false;
        }

        bool empty() const { return m_regions.empty(); }

        std::set< Region >::const_iterator cbegin() const { return m_regions.cbegin(); }
        std::set< Region >::const_iterator cend() const { return m_regions.cend(); }

        std::set< Region >::iterator begin() const { return m_regions.begin(); }
        std::set< Region >::iterator end() const { return m_regions.end(); }

        std::string toString() const;
        // Utility function so we can std::cout a region
        friend std::ostream & operator<<( std::ostream & out, const SetRegions & theRegions );

        std::size_t size() const { return m_regions.size(); }

    private:
        std::set< Region > m_regions;
    };

    using regions_t = std::vector< caller::Region >;
}
}

#endif
