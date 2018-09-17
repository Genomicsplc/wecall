// All content Copyright (C) 2018 Genomics plc
//
// Created by adorr on 01/09/17.
//

#ifndef WECALL_MM_HELPERS_H
#define WECALL_MM_HELPERS_H

#include <emmintrin.h>
#include <cstddef>

struct short_array8
{
public:
    short_array8( short value ) : m_value( _mm_set1_epi16( value ) ) {}
    short_array8( short v0, short v1, short v2, short v3, short v4, short v5, short v6, short v7 )
        : m_value( _mm_set_epi16( v7, v6, v5, v4, v3, v2, v1, v0 ) )
    {
    }
    short_array8() : m_value() {}
    explicit short_array8( __m128i value ) : m_value( value ) {}
    short_array8( const short_array8 & value ) : m_value( value.m_value ) {}
    short operator[]( std::size_t index ) const { return ( (const int16_t *)&m_value )[index]; }
    template < typename FUNC >
    short_array8( FUNC f )
        : short_array8( f( 0 ), f( 1 ), f( 2 ), f( 3 ), f( 4 ), f( 5 ), f( 6 ), f( 7 ) )
    {
    }

    class reference
    {
    public:
        reference( short_array8 & array, short index ) : m_array( array ), m_index( index ) {}
        reference & operator=( short value )
        {
            m_array.m_value = _mm_insert_epi16( m_array.m_value, value, m_index );
            return *this;
        }
        operator short() const
        {
            const auto & array = m_array;
            return array[m_index];
        }

    private:
        short_array8 & m_array;
        short m_index;
    };
    reference operator[]( std::size_t index ) { return reference( *this, index ); }
    short_array8 & operator=( const short_array8 & other )
    {
        m_value = other.m_value;
        return *this;
    }

    __m128i m_value;
};

short_array8 min( const short_array8 & v1, const short_array8 & v2 )
{
    return short_array8( _mm_min_epi16( v1.m_value, v2.m_value ) );
}

short_array8 andnot( const short_array8 & v1, const short_array8 & v2 )
{
    return short_array8( _mm_andnot_si128( v1.m_value, v2.m_value ) );
}

short_array8 cmpeq( const short_array8 & v1, const short_array8 & v2 )
{
    return short_array8( _mm_cmpeq_epi16( v1.m_value, v2.m_value ) );
}

short_array8 operator+( const short_array8 & v1, const short_array8 & v2 )
{
    return short_array8( _mm_add_epi16( v1.m_value, v2.m_value ) );
}

short_array8 operator&( const short_array8 & v1, const short_array8 & v2 )
{
    return short_array8( _mm_and_si128( v1.m_value, v2.m_value ) );
}

short_array8 operator|( const short_array8 & v1, const short_array8 & v2 )
{
    return short_array8( _mm_or_si128( v1.m_value, v2.m_value ) );
}

short_array8 operator<<( const short_array8 & v1, int numbits )
{
    return short_array8( _mm_slli_epi16( v1.m_value, numbits ) );
}

short_array8 operator>>( const short_array8 & v1, int numbits )
{
    return short_array8( _mm_srli_epi16( v1.m_value, numbits ) );
}

short_array8 shift_right( const short_array8 & v ) { return short_array8( _mm_srli_si128( v.m_value, 2 ) ); }

short_array8 shift_left( const short_array8 & v ) { return short_array8( _mm_slli_si128( v.m_value, 2 ) ); }

#endif  // WECALL_MM_HELPERS_H
