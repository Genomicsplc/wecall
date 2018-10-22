// All content Copyright (C) 2018 Genomics plc
#include <sstream>
#include <vector>

#include "utils/interval.hpp"

namespace wecall
{
namespace utils
{
    bool Interval::overlapsOrTouches( const Interval & other ) const
    {
        return this->overlaps( other ) or this->end() == other.start() or this->start() == other.end();
    }

    Interval Interval::getIntersect( const Interval & other ) const
    {
        WECALL_ASSERT( overlaps( other ),
                        "Cant intersect non-overlapping intervals: " + toString() + " & " + other.toString() );
        return Interval( std::max( m_start, other.m_start ), std::min( m_end, other.m_end ) );
    }

    void Interval::combine( const Interval & other )
    {
        m_start = std::min( m_start, other.m_start );
        m_end = std::max( m_end, other.m_end );
    }

    Interval Interval::getPadded( const int64_t & paddingSize ) const
    {
        return Interval( m_start - paddingSize, m_end + paddingSize );
    }

    Interval & Interval::operator+=( const int64_t offset )
    {
        this->m_start += offset;
        this->m_end += offset;
        return *this;
    }
    Interval & Interval::operator-=( const int64_t offset )
    {
        this->m_start -= offset;
        this->m_end -= offset;
        return *this;
    }

    std::string Interval::toString() const
    {
        std::stringstream sstrInterval;
        sstrInterval << *this;
        return sstrInterval.str();
    }

    std::ostream & operator<<( std::ostream & out, const Interval & theInterval )
    {
        return out << "[" << theInterval.start() << ", " << theInterval.end() << ")";
    }

    bool Interval::operator<( const Interval & other ) const
    {
        if ( m_start == other.m_start )
        {
            return m_end < other.m_end;
        }

        return m_start < other.m_start;
    }
}
}
