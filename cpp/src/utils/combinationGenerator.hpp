// All content Copyright (C) 2018 Genomics plc
#ifndef _COMBINATION_GENERATOR_H_
#define _COMBINATION_GENERATOR_H_
#include <vector>

namespace wecall
{
class CombinationIterator
{
public:
    CombinationIterator( const long int n, const long int m )
        : m_n( n ), m_m( m ), m_indices( n, 0u ), m_carry( m == 0 )
    {
    }

    CombinationIterator() : m_n( 0 ), m_m( 0 ), m_indices(), m_carry( true ) {}

    const std::vector< unsigned long > & operator*() const { return m_indices; }

    CombinationIterator & operator++()
    {
        while ( not m_carry )
        {
            ++m_indices.back();
            for ( auto it = m_indices.rbegin(); it != m_indices.rend(); ++it )
            {
                if ( m_carry )
                {
                    ( *it )++;
                    m_carry = false;
                }
                if ( *it == m_m )
                {
                    *it = 0;
                    m_carry = true;
                }
                if ( not m_carry )
                {
                    break;
                }
            }
            unsigned long max_elem = 0;
            for ( const auto & index : m_indices )
            {
                if ( index < max_elem )
                {
                    max_elem = m_m;
                    break;
                }
                max_elem = std::max( max_elem, index );
            }
            if ( max_elem < m_m )
            {
                break;
            }
        }
        return *this;
    }

    bool operator!=( const CombinationIterator & rhs ) const { return not m_carry; }

private:
    const unsigned long m_n;
    const unsigned long m_m;
    std::vector< unsigned long > m_indices;
    bool m_carry = false;
};

class CombinationGenerator
{
public:
    CombinationGenerator( const unsigned long n, const unsigned long m ) : m_n( n ), m_m( m ) {}

    CombinationIterator begin() const { return CombinationIterator( m_n, m_m ); }
    CombinationIterator end() const { return CombinationIterator(); }

private:
    const unsigned long m_n;
    const unsigned long m_m;
};
}

#endif
