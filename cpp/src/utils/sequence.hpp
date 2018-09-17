// All content Copyright (C) 2018 Genomics plc
#ifndef WECALL_SEQUENCE_HPP_H
#define WECALL_SEQUENCE_HPP_H

#include <string>
#include <memory>

#include "boost/operators.hpp"

#include "common.hpp"

namespace echidna
{
namespace utils
{
    class BasePairSequence
        : boost::partially_ordered< BasePairSequence, BasePairSequence, boost::equality_comparable< BasePairSequence > >
    {
    public:
        BasePairSequence() : m_sequence() {}
        BasePairSequence( const char * sequence ) : m_sequence( sequence ) {}
        BasePairSequence( const std::string & sequence ) : m_sequence( sequence ) {}
        BasePairSequence( std::size_t n, char c ) : m_sequence( std::string( n, c ) ) {}
        BasePairSequence( const BasePairSequence & other, std::size_t pos, std::size_t length );
        BasePairSequence( const BasePairSequence & other ) : m_sequence( other.m_sequence ) {}
        bool operator==( const BasePairSequence & other ) const;
        bool operator<( const BasePairSequence & other ) const;
        BasePairSequence & operator=( const BasePairSequence & other )
        {
            m_sequence = other.m_sequence;
            return *this;
        }

        char operator[]( std::size_t pos ) const { return at( pos ); }
        const char at( std::size_t pos ) const { return m_sequence.at( pos ); }

        std::string str() const { return m_sequence; };
        std::size_t size() const { return m_sequence.size(); }

        BasePairSequence substr( int64_t pos, int64_t length ) const;
        BasePairSequence substr( int64_t pos ) const;
        char front() const { return m_sequence.front(); }
        char back() const { return m_sequence.back(); }

        bool contains( const char c ) const;

        using const_iterator = std::string::const_iterator;
        using const_reverse_iterator = std::string::const_reverse_iterator;

        const_iterator cbegin() const { return m_sequence.cbegin(); }
        const_iterator cend() const { return m_sequence.cend(); }
        const_reverse_iterator crbegin() const { return m_sequence.crbegin(); }
        const_reverse_iterator crend() const { return m_sequence.crend(); }

        friend std::ostream & operator<<( std::ostream & out, const BasePairSequence & __str )
        {
            return out << __str.str();
        }

        BasePairSequence operator+( const BasePairSequence & other ) const;

        BasePairSequence leftTrimmed() const;
        BasePairSequence rightTrimmed() const;

    private:
        std::string m_sequence;
    };

    using QualitySequence = std::string;
}
}

namespace std
{
template <>
struct hash< echidna::utils::BasePairSequence >
{
    std::size_t operator()( const echidna::utils::BasePairSequence & x ) const
    {
        return hash< std::string >()( x.str() );
    }
};
}

#endif  // WECALL_SEQUENCE_HPP_H
