// All content Copyright (C) 2018 Genomics plc
#include "utils/referenceSequence.hpp"
#include "caller/region.hpp"
#include "utils/sequence.hpp"

namespace echidna
{
namespace utils
{
    ReferenceSequence::ReferenceSequence( caller::Region region, echidna::utils::BasePairSequence sequence )
        : m_region( region ), m_sequence( sequence )
    {
        ECHIDNA_ASSERT( ( m_region.size() == m_sequence.size() ),
                        "ReferenceSequence requires region size to match sequence size: " +
                            std::to_string( m_region.size() ) + " != " + std::to_string( m_sequence.size() ) );
    }

    ReferenceSequence ReferenceSequence::subseq( const caller::Region & region ) const
    {
        ECHIDNA_ASSERT( m_region.contains( region ),
                        "Cannot get subsequence for region not contained in reference sequence: " + region.toString() +
                            " must be contained " + m_region.toString() );
        return ReferenceSequence( region, m_sequence.substr( distanceFromStart( region.start() ), region.size() ) );
    }

    std::pair< std::string::const_iterator, std::string::const_iterator > ReferenceSequence::subseqRange(
        const caller::Region & region ) const
    {
        ECHIDNA_ASSERT( m_region.contains( region ),
                        "Cannot get subsequence for region not contained in reference sequence: " + region.toString() +
                            " must be contained " + m_region.toString() );

        auto it = m_sequence.cbegin() + distanceFromStart( region.start() );
        return std::make_pair( it, it + region.size() );
    }

    std::pair< ReferenceSequence::const_iterator, ReferenceSequence::const_iterator >
    ReferenceSequence::getRangeForwardIterators( utils::Interval interval ) const
    {
        ECHIDNA_ASSERT( m_region.interval().contains( interval ),
                        "Interval " + interval.toString() + " not contained in " + m_region.toString() );
        const auto start = m_sequence.cbegin() + int64_to_sizet( distanceFromStart( interval.start() ) );
        const auto end = start + int64_to_sizet( interval.size() );
        return std::make_pair( start, end );
    }

    std::pair< ReferenceSequence::const_reverse_iterator, ReferenceSequence::const_reverse_iterator >
    ReferenceSequence::getRangeReverseIterators( utils::Interval interval ) const
    {
        ECHIDNA_ASSERT( m_region.interval().contains( interval ),
                        "Interval " + interval.toString() + " not contained in " + m_region.toString() );
        const auto start = m_sequence.crbegin() + int64_to_sizet( distanceFromEnd( interval.end() ) );
        const auto end = start + int64_to_sizet( interval.size() );
        return std::make_pair( start, end );
    }

    ReferenceSequence ReferenceSequence::getPadded( const caller::Region & widerRegion ) const
    {
        ECHIDNA_ASSERT( widerRegion.contains( m_region ),
                        "Cannot get padded sequence for region not containing original region: " +
                            widerRegion.toString() + " must contain " + m_region.toString() );

        const auto newLengthAtStart = m_region.start() - widerRegion.start();
        const auto newLengthAtEnd = widerRegion.end() - m_region.end();
        const BasePairSequence leftPadding( int64_to_sizet( newLengthAtStart ), constants::gapChar );
        const BasePairSequence rightPadding( int64_to_sizet( newLengthAtEnd ), constants::gapChar );

        return ReferenceSequence( widerRegion, leftPadding + m_sequence + rightPadding );
    }

    ReferenceSequence ReferenceSequence::getPaddedSequence( const caller::Region & clusterRegion,
                                                            const caller::Region & maximumRefFileRegion,
                                                            const int64_t desiredPadding ) const
    {
        ECHIDNA_ASSERT( this->region().contains( maximumRefFileRegion ),
                        "Maximal padding region must be contained in reference." );
        ECHIDNA_ASSERT( maximumRefFileRegion.contains( clusterRegion ),
                        "Maximal cluster region " + maximumRefFileRegion.toString() +
                            " does not contain cluster region " + clusterRegion.toString() );

        const auto desiredRegion = clusterRegion.getPadded( desiredPadding );
        const auto refFileRegion = desiredRegion.getIntersect( maximumRefFileRegion );
        const auto reference = this->subseq( refFileRegion );
        return reference.getPadded( desiredRegion );
    }
}
}
