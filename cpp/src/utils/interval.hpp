// All content Copyright (C) 2018 Genomics plc
#ifndef INTERVAL_HPP
#define INTERVAL_HPP

#include <string>
#include <vector>
#include "utils/logging.hpp"

namespace wecall
{
namespace utils
{
    /// struct interval to handle intervals internally. Logs and throws when non-proper intervals are attempted to be
    /// created.
    class Interval
    {
    public:
        /// Construct an interval from left and right bounds.
        Interval( int64_t _left, int64_t _right ) : m_start( _left ), m_end( _right )
        {
            ECHIDNA_ASSERT( m_start <= m_end, this->toString() + " is not valid (half-open) interval" );
        }

        /// @return The centre of the interval. Takes floor in integral types.
        int64_t centre() const { return ( m_start + m_end ) / 2L; }

        int64_t size() const { return m_end - m_start; }

        bool contains( int64_t otherPos ) const { return m_start <= otherPos and otherPos < m_end; }
        bool contains( const Interval & other ) const { return other.m_start >= m_start and other.m_end <= m_end; }
        bool overlaps( const Interval & other ) const
        {
            if ( *this == other )
                return true;

            return m_start < other.m_end and m_end > other.m_start;
        }

        bool overlapsOrTouches( const Interval & other ) const;

        void combine( const Interval & other );

        Interval getPadded( const int64_t & paddingSize ) const;
        Interval getIntersect( const Interval & other ) const;

        Interval & operator+=( const int64_t offset );
        Interval & operator-=( const int64_t offset );

        Interval operator+( const int64_t offset ) const { return Interval( *this ) += offset; }
        Interval operator-( const int64_t offset ) const { return Interval( *this ) -= offset; }

        bool operator==( const Interval & other ) const { return m_start == other.m_start and m_end == other.m_end; }

        /// @return True if other is larger using a lexigraphical ordering.
        bool operator<( const Interval & other ) const;

        std::string toString() const;

        // Utility function so we can std::cout a region
        friend std::ostream & operator<<( std::ostream & out, const Interval & theInterval );

        int64_t start() const { return m_start; }
        int64_t end() const { return m_end; }

    private:
        int64_t m_start;
        int64_t m_end;
    };
}
}
#endif
