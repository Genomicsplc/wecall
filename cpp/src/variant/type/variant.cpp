// All content Copyright (C) 2018 Genomics plc
#include "common.hpp"
#include "caller/region.hpp"
#include "utils/interval.hpp"
#include "variant/type/variant.hpp"

#include <algorithm>
#include <sstream>
#include <boost/optional.hpp>

namespace echidna
{
namespace variant
{

    Variant::Variant( utils::referenceSequencePtr_t referenceSequencePtr,
                      Region region,
                      BasePairSequence added,
                      bool isFullyLeftAligned )
        : m_refSequence( referenceSequencePtr ),
          m_region( region ),
          m_altSequence( added ),
          m_isFullyLeftAligned( isFullyLeftAligned ),
          m_prior( std::numeric_limits< double >::quiet_NaN() ),
          m_neverFilter( false ),
          m_isGenotypingVariant( false ),
          m_fromBreakpoint( false )
    {
    }

    bool Variant::operator==( const Variant & other ) const
    {
        return refSequence() == other.refSequence() and sequence() == other.sequence();
    }

    caller::Region Variant::getLeftAlignedRegion( const boost::optional< int64_t > minPos ) const
    {
        const int64_t nRemoved = m_region.size();

        if ( not this->isPureIndel() )
        {
            return this->region();
        }

        auto minPosition = minPos.is_initialized() ? minPos.get() : m_refSequence->region().interval().start();

        const auto refIterators =
            m_refSequence->getRangeReverseIterators( utils::Interval( minPosition, m_region.end() ) );

        auto refIt = refIterators.first;

        for ( auto mutIt( m_altSequence.crbegin() ); refIt != refIterators.second and mutIt != m_altSequence.crend();
              ++mutIt, ++refIt )
        {
            if ( ( *refIt ) != ( *mutIt ) )
            {
                const auto newRegionStart = minPosition + std::distance( refIt, refIterators.second ) - nRemoved;
                return caller::Region( this->contig(), newRegionStart, newRegionStart + nRemoved );
            }
        }

        const auto refBeforeIterators =
            m_refSequence->getRangeReverseIterators( utils::Interval( minPosition, m_region.start() ) );

        for ( auto mutIt( refBeforeIterators.first );
              refIt != refIterators.second and mutIt != refBeforeIterators.second; ++mutIt, ++refIt )
        {
            if ( ( *refIt ) != ( *mutIt ) )
            {
                const auto newRegionStart = minPosition + std::distance( refIt, refIterators.second ) - nRemoved;
                return caller::Region( this->contig(), newRegionStart, newRegionStart + nRemoved );
            }
        }

        ECHIDNA_LOG( SUPER_DEBUG, "Could not completely left align variant: " << this->toString() << " on reference "
                                                                              << m_refSequence->region()
                                                                              << " with minPos " << minPosition );

        return caller::Region( m_region.contig(), minPosition, minPosition + nRemoved );
    }

    caller::Region Variant::getRightAlignedRegion( const boost::optional< int64_t > maxPos ) const
    {
        const int64_t nRemoved = m_region.size();

        if ( not this->isPureIndel() )
        {
            return this->region();
        }

        auto maxPosition = maxPos.is_initialized() ? maxPos.get() : m_refSequence->region().interval().end();

        const auto refIterators =
            m_refSequence->getRangeForwardIterators( utils::Interval( m_region.start(), maxPosition ) );

        auto refIt = refIterators.first;

        for ( auto mutIt( m_altSequence.cbegin() ); refIt != refIterators.second and mutIt != m_altSequence.cend();
              ++mutIt, ++refIt )
        {
            if ( ( *refIt ) != ( *mutIt ) )
            {
                const auto newRegionStart = m_region.start() + std::distance( refIterators.first, refIt );
                return caller::Region( this->contig(), newRegionStart, newRegionStart + nRemoved );
            }
        }

        const auto refAfterIterators =
            m_refSequence->getRangeForwardIterators( utils::Interval( m_region.end(), maxPosition ) );

        for ( auto mutIt( refAfterIterators.first ); refIt != refIterators.second and mutIt != refAfterIterators.second;
              ++mutIt, ++refIt )
        {
            if ( ( *refIt ) != ( *mutIt ) )
            {
                const auto newRegionStart = m_region.start() + std::distance( refIterators.first, refIt );
                return caller::Region( this->contig(), newRegionStart, newRegionStart + nRemoved );
            }
        }

        ECHIDNA_LOG( SUPER_DEBUG, "Could not completely right align variant: " << this->toString() << " on reference "
                                                                               << m_refSequence->region()
                                                                               << " with maxPos " << maxPosition );
        return caller::Region( m_region.contig(), maxPosition - nRemoved, maxPosition );
    }

    varPtr_t Variant::getLeftAligned( const int64_t minPos ) const
    {
        const int64_t nAdded = static_cast< int64_t >( m_altSequence.size() );

        const auto newRegion = this->getLeftAlignedRegion( minPos );
        const auto fullyLeftAligned = ( newRegion.start() != minPos );

        if ( newRegion == m_region )
        {
            return std::make_shared< Variant >( *this );
        }
        else if ( newRegion.start() + nAdded <= m_region.start() )
        {
            const auto mutSequence = m_refSequence->subseq( caller::Region( this->contig(), newRegion.start(),
                                                                            newRegion.start() + nAdded ) ).sequence();
            return std::make_shared< Variant >( m_refSequence, newRegion, mutSequence, fullyLeftAligned );
        }
        else
        {
            const auto mutSequence = m_refSequence->subseq( caller::Region( this->contig(), newRegion.start(),
                                                                            m_region.start() ) ).sequence() +
                                     m_altSequence.substr( 0, newRegion.start() + nAdded - m_region.start() );
            return std::make_shared< Variant >( m_refSequence, newRegion, mutSequence, fullyLeftAligned );
        }
    }

    varPtr_t Variant::getRightAligned( const int64_t maxPos ) const
    {
        const int64_t nAdded = m_altSequence.size();

        const auto newRegion = this->getRightAlignedRegion( maxPos );

        if ( newRegion == m_region )
        {
            return std::make_shared< Variant >( *this );
        }
        else if ( newRegion.end() - nAdded >= m_region.end() )
        {
            const auto mutSequence = m_refSequence->subseq( caller::Region( this->contig(), newRegion.end() - nAdded,
                                                                            newRegion.end() ) ).sequence();
            return std::make_shared< Variant >( m_refSequence, newRegion, mutSequence, false );
        }
        else
        {
            const auto length = nAdded - ( newRegion.end() - m_region.end() );
            const auto mutSequence =
                m_altSequence.substr( newRegion.end() - m_region.end(), length ) +
                m_refSequence->subseq( caller::Region( this->contig(), m_region.end(), newRegion.end() ) ).sequence();
            return std::make_shared< Variant >( m_refSequence, newRegion, mutSequence, false );
        }
    }

    varPtr_t Variant::getTrimmed() const
    {
        const auto refSequence = m_refSequence->subseq( m_region );
        size_t tailLength = 0;
        for ( auto refIt( refSequence.crbegin() ), mutIt( m_altSequence.crbegin() );
              refIt != refSequence.crend() and mutIt != m_altSequence.crend() and *refIt == *mutIt; ++mutIt, ++refIt )
        {
            ++tailLength;
        }

        size_t headLength = 0;
        for ( auto refIt( refSequence.cbegin() ), mutIt( m_altSequence.cbegin() );
              refIt != refSequence.cend() and mutIt != m_altSequence.cend() and *refIt == *mutIt; ++mutIt, ++refIt )
        {
            ++headLength;
        }

        const auto totalTrim = static_cast< int >( tailLength + headLength );

        const auto removedOverlap = totalTrim - static_cast< int >( m_region.size() );
        const auto addedOverlap = totalTrim - static_cast< int >( m_altSequence.size() );

        size_t adjustedHeadLength = headLength;
        if ( removedOverlap > 0 or addedOverlap > 0 )
        {
            adjustedHeadLength = headLength - static_cast< size_t >( std::max( removedOverlap, addedOverlap ) );
        }

        caller::Region region( this->contig(), this->start() + adjustedHeadLength, this->end() - tailLength );
        const auto newSeq =
            m_altSequence.substr( adjustedHeadLength, m_altSequence.size() - adjustedHeadLength - tailLength );

        return std::make_shared< Variant >( m_refSequence, region, newSeq );
    }

    Variant Variant::operator-( const Variant & other ) const
    {
        ECHIDNA_ASSERT( this->region().contains( other.region() ), "" );
        ECHIDNA_ASSERT( this->m_altSequence.size() >= other.m_altSequence.size(), "" );
        if ( other.end() == this->end() )
        {
            const auto newSeq = m_altSequence.substr( 0, m_altSequence.size() - other.m_altSequence.size() );
            return Variant( m_refSequence, Region( this->contig(), this->start(), other.start() ), newSeq );
        }
        else if ( other.start() == this->start() )
        {
            const auto newSeq =
                m_altSequence.substr( other.m_altSequence.size(), m_altSequence.size() - other.m_altSequence.size() );
            return Variant( m_refSequence, Region( this->contig(), other.end(), this->end() ), newSeq );
        }
        else
        {
            throw utils::echidna_exception( "Only can substract variants starting or ending at same place" );
        }
    }

    Variant Variant::operator+( const Variant & other ) const
    {
        ECHIDNA_ASSERT( m_refSequence->region().contains( other.region() ), "" );
        ECHIDNA_ASSERT( m_region.end() == other.start(), "Only can join neighbours" );
        Region newRegion = m_region;
        newRegion.combine( other.region() );
        return Variant( m_refSequence, newRegion, m_altSequence + other.m_altSequence );
    }

    varPtr_t Variant::remove( const varPtr_t & other ) const
    {
        ECHIDNA_ASSERT( this->removable( other ), "Only can remove removable variants" );
        return ( *this - *other ).getTrimmed();
    }

    bool Variant::removable( const varPtr_t & other ) const
    {
        if ( not this->region().contains( other->region() ) )
        {
            return false;
        }
        else if ( this->m_altSequence.size() < other->m_altSequence.size() )
        {
            return false;
        }
        else if ( this->end() == other->end() )
        {
            // Remove from end
            const auto altEnd = this->m_altSequence.substr( this->m_altSequence.size() - other->m_altSequence.size(),
                                                            other->m_altSequence.size() );
            return altEnd == other->m_altSequence;
        }
        else if ( this->start() == other->start() )
        {
            // Remove at start
            const auto altStart = this->m_altSequence.substr( 0, other->m_altSequence.size() );
            return altStart == other->m_altSequence;
        }
        else
        {
            // Currently not removable
            return false;
        }
    }

    bool Variant::joinable( const varPtr_t & other ) const
    {
        return this->contig() == other->contig() and this->end() == other->start();
    }

    std::vector< varPtr_t > Variant::split() const
    {
        if ( empty() )
        {
            return {};
        }
        else if ( isSNP() or isPureIndel() )
        {
            return {std::make_shared< Variant >( *this )};
        }
        else if ( sequenceLength() == sequenceLengthInRef() )
        {
            std::vector< varPtr_t > splitMNP;
            const auto referenceSequence = this->refSequence().sequence();
            for ( std::size_t i = 0; i < static_cast< size_t >( sequenceLengthInRef() ); ++i )
            {
                if ( m_altSequence[i] != referenceSequence[i] )
                {
                    const Region region( contig(), start() + i, start() + i + 1 );
                    splitMNP.push_back(
                        std::make_shared< Variant >( m_refSequence, region, m_altSequence.substr( i, 1 ) ) );
                }
            }
            return splitMNP;
        }
        else if ( isDeletion() )
        {
            const auto pureDeletionEnd = end() - sequence().size();

            std::vector< varPtr_t > split = {
                std::make_shared< Variant >( m_refSequence, Region( contig(), start(), pureDeletionEnd ), "" )};

            const auto mnpRegion = Region( contig(), pureDeletionEnd, end() );
            assert( std::size_t( mnpRegion.size() ) == m_altSequence.size() );
            const auto splitMNPRemaining = Variant( m_refSequence, mnpRegion, m_altSequence ).split();
            split.insert( split.end(), splitMNPRemaining.cbegin(), splitMNPRemaining.cend() );

            return split;
        }
        else if ( isInsertion() )
        {
            const auto altIns = m_altSequence.substr( 0, sequenceLengthChange() );
            const auto mnpAdded =
                m_altSequence.substr( sequenceLengthChange(), m_altSequence.size() - sequenceLengthChange() );

            std::vector< varPtr_t > split = {
                std::make_shared< Variant >( m_refSequence, Region( contig(), start(), start() ), altIns )};

            assert( std::size_t( m_region.size() ) == mnpAdded.size() );
            const auto splitMNPRemaining = Variant( m_refSequence, m_region, mnpAdded ).split();
            split.insert( split.end(), splitMNPRemaining.cbegin(), splitMNPRemaining.cend() );

            return split;
        }
        else
        {
            return {std::make_shared< Variant >( *this )};
        }
    }

    void setDefaultPriors( const std::vector< varPtr_t > & variants )
    {
        for ( const auto & var : variants )
        {
            if ( std::isnan( var->prior() ) )
            {
                if ( var->isDeletion() )
                {
                    var->prior( 1e-4 * pow( 0.8, var->sequenceLengthInRef() ) );
                }
                else if ( var->isInsertion() )
                {
                    var->prior( 1e-4 * pow( 0.33, var->sequenceLength() ) );
                }
                else if ( var->isSNP() )
                {
                    var->prior( 1.0e-3 / 3.0 );
                }
                else if ( var->sequenceLength() == var->sequenceLengthInRef() )
                {
                    // Gerton Lunter's setPrior for MNPs and length-preserving replacements.
                    auto nDiffs = 0;

                    for ( std::size_t i = 0; i < static_cast< size_t >( var->sequenceLengthInRef() ); ++i )
                    {
                        if ( var->sequence()[i] != var->refSequence().sequence()[i] )
                        {
                            ++nDiffs;
                        }
                    }

                    var->prior( 5e-5 * pow( 0.1, nDiffs - 1 ) * ( 1.0 - 0.1 ) );
                }
                else
                {
                    //
                }
            }
        }
    }

    std::string Variant::toString() const
    {
        std::stringstream repr;
        repr << "Variant(" << this->refSequence() << " --> " << m_altSequence << ")";
        return repr.str();
    }

    caller::SetRegions Variant::getStartEndRegions( const boost::optional< int64_t > minPos,
                                                    const boost::optional< int64_t > maxPos ) const
    {
        const auto leftAligned = getLeftAlignedRegion( minPos );
        const auto rightAligned = getRightAlignedRegion( maxPos );
        const auto leftRegion = caller::Region( this->contig(), leftAligned.start(), rightAligned.start() );
        const auto rightRegion = caller::Region( this->contig(), leftAligned.end(), rightAligned.end() );

        caller::SetRegions setRegions;
        setRegions.insert( leftRegion );
        setRegions.insert( rightRegion );
        return setRegions;
    }

    //-------------------------------------------------------------------------------------------------

    bool varPtrComp::operator()( const varPtr_t & x, const varPtr_t & y ) const
    {
        if ( x->contig() != y->contig() )
        {
            return x->contig() < y->contig();
        }

        if ( x->start() != y->start() )
        {
            return x->start() < y->start();
        }

        if ( x->zeroIndexedVcfPosition() != y->zeroIndexedVcfPosition() )
        {
            return x->zeroIndexedVcfPosition() < y->zeroIndexedVcfPosition();
        }

        if ( x->end() != y->end() )
        {
            return x->end() < y->end();
        }

        if ( x->sequenceLength() != y->sequenceLength() )
        {
            return x->sequenceLength() < y->sequenceLength();
        }

        return x->sequence() < y->sequence();
    }
}
}
