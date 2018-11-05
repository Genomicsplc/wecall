// All content Copyright (C) 2018 Genomics plc
#include "utils/sequence.hpp"

namespace wecall
{
namespace utils
{
    template < typename Iterator >
    std::size_t findFirstNormalChar( Iterator begin, Iterator end )
    {
        std::size_t index = 0;
        while ( begin != end and *begin == constants::gapChar )
        {
            ++begin;
            ++index;
        }
        return index;
    }

    BasePairSequence BasePairSequence::substr( int64_t pos, int64_t length ) const
    {
        if ( pos < 0 )
            throw wecall::utils::wecall_exception( "BasePairSequence: pos less than 0" );
        return BasePairSequence( m_sequence.substr( pos, length ) );
    }

    BasePairSequence BasePairSequence::substr( int64_t pos ) const
    {
        if ( pos < 0 )
            throw wecall::utils::wecall_exception( "BasePairSequence: pos less than 0" );
        return BasePairSequence( m_sequence.substr( pos ) );
    }

    BasePairSequence BasePairSequence::operator+( const BasePairSequence & other ) const
    {
        return BasePairSequence( str() + other.str() );
    }

    BasePairSequence BasePairSequence::leftTrimmed() const
    {
        const std::size_t index = findFirstNormalChar( this->cbegin(), this->cend() );
        return this->substr( index, this->size() - index );
    }

    BasePairSequence BasePairSequence::rightTrimmed() const
    {
        const std::size_t index = findFirstNormalChar( this->crbegin(), this->crend() );
        return this->substr( 0, this->size() - index );
    }

    bool BasePairSequence::operator==( const BasePairSequence & other ) const { return m_sequence == other.m_sequence; }

    bool BasePairSequence::operator<( const BasePairSequence & other ) const { return m_sequence < other.m_sequence; }

    bool BasePairSequence::contains( const char c ) const { return m_sequence.find( c ) != std::string::npos; }
}
}
