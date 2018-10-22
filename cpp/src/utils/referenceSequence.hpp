// All content Copyright (C) 2018 Genomics plc
#ifndef ALIGNED_SEQUENCE_HPP
#define ALIGNED_SEQUENCE_HPP

#include <string>

#include "common.hpp"
#include "utils/logging.hpp"
#include "utils/sequence.hpp"
#include "caller/region.hpp"

namespace wecall
{
namespace utils
{
    class ReferenceSequence
    {
    public:
        ReferenceSequence( caller::Region region, utils::BasePairSequence sequence );

        bool operator==( const ReferenceSequence & other ) const
        {
            return m_region == other.m_region and m_sequence == other.m_sequence;
        }

        const caller::Region & region() const { return m_region; }
        const utils::BasePairSequence & sequence() const { return m_sequence; }

        char at( int64_t refPos ) const { return m_sequence.at( int64_to_sizet( distanceFromStart( refPos ) ) ); }

        std::size_t size() const { return m_sequence.size(); }
        std::string contig() const { return m_region.contig(); }
        int64_t start() const { return m_region.start(); }
        int64_t end() const { return m_region.end(); }

        ReferenceSequence subseq( const caller::Region & region ) const;
        std::pair< std::string::const_iterator, std::string::const_iterator > subseqRange(
            const caller::Region & region ) const;

        ReferenceSequence getPadded( const caller::Region & widerRegion ) const;

        // Return reference sequence with pads outside the maximal region with 'N's but returns sequence padded by at
        // least the desired padding.
        ReferenceSequence getPaddedSequence( const caller::Region & clusterRegion,
                                             const caller::Region & maximumRefFileRegion,
                                             const int64_t desiredPadding ) const;

        friend std::ostream & operator<<( std::ostream & out, const ReferenceSequence & referenceSequence )
        {
            return out << referenceSequence.region() << " " << referenceSequence.sequence().str();
        }

        std::string toString() const { return m_region.toString() + "\t" + m_sequence.str(); }

        using const_iterator = wecall::utils::BasePairSequence::const_iterator;
        using const_reverse_iterator = wecall::utils::BasePairSequence::const_reverse_iterator;

        const_iterator cbegin() const { return m_sequence.cbegin(); }
        const_iterator cend() const { return m_sequence.cend(); }

        const_reverse_iterator crbegin() const { return m_sequence.crbegin(); }
        const_reverse_iterator crend() const { return m_sequence.crend(); }

        std::pair< const_iterator, const_iterator > getRangeForwardIterators( utils::Interval interval ) const;
        std::pair< const_reverse_iterator, const_reverse_iterator > getRangeReverseIterators(
            utils::Interval interval ) const;

    private:
        int64_t distanceFromStart( int64_t refPos ) const { return refPos - m_region.start(); }
        int64_t distanceFromEnd( int64_t refPos ) const { return m_region.end() - refPos; }

        caller::Region m_region;
        utils::BasePairSequence m_sequence;
    };

    typedef std::shared_ptr< ReferenceSequence > referenceSequencePtr_t;
}
}

#endif
