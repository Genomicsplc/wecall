// All content Copyright (C) 2018 Genomics plc
#include "samtools/bam.h"
#include <utils/exceptions.hpp>
#include <algorithm>
#include <iostream>
#include "alignment/cigar.hpp"
#include "alignment/cigarItems.hpp"

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <unordered_map>
#include <numeric>

namespace wecall
{
namespace alignment
{
    static int64_t cigarItemCacheSize = 250;
    static std::shared_ptr< CigarItem > rollCigarItem( cigarFlags flag, int64_t length )
    {
        switch ( flag )
        {

        case cigarFlags::MATCH:
            return std::make_shared< CigarMatch >( length );
        case cigarFlags::INSERTION:
            return std::make_shared< CigarInsertion >( length );
        case cigarFlags::DELETION:
            return std::make_shared< CigarDeletion >( length );
        case cigarFlags::SKIP:
            return std::make_shared< CigarSkip >( length );
        case cigarFlags::SOFT_CLIP:
            return std::make_shared< CigarSoftClip >( length );
        case cigarFlags::HARD_CLIP:
            return std::make_shared< CigarHardClip >( length );
        case cigarFlags::PADDING:
            return std::make_shared< CigarPadding >( length );
        case cigarFlags::SEQ_MATCH:
            return std::make_shared< CigarSequenceMatch >( length );
        case cigarFlags::SEQ_MISMATCH:
            return std::make_shared< CigarSequenceMismatch >( length );
        default:
            throw utils::wecall_exception( "Cigar flag number: " + std::to_string( static_cast< int >( flag ) ) +
                                            " not implemented." );
        }
    }

    struct pairhash
    {
    public:
        std::size_t operator()( const std::pair< cigarFlags, int64_t > & x ) const
        {
            return std::hash< int >()( int( x.first ) ) ^ std::hash< int64_t >()( x.second );
        }
    };

    typedef std::unordered_map< std::pair< cigarFlags, int64_t >, std::shared_ptr< CigarItem >, pairhash >
        cigarItemMap_t;

    static cigarItemMap_t makeCigarItemCache()
    {
        cigarItemMap_t ret;
        for ( int64_t i = 0; i < cigarItemCacheSize; ++i )
        {
            for ( auto j = int( cigarFlags::MATCH ); j <= int( cigarFlags::SEQ_MISMATCH ); j++ )
            {
                auto flag = cigarFlags( j );
                ret[std::make_pair( flag, i )] = rollCigarItem( flag, i );
            }
        }
        return ret;
    }

    static auto cigarItemCache = makeCigarItemCache();

    std::shared_ptr< CigarItem > CigarItemBuilder::roll( cigarFlags flag, int64_t length )
    {
        switch ( flag )
        {

        case cigarFlags::MATCH:
        case cigarFlags::INSERTION:
        case cigarFlags::DELETION:
        case cigarFlags::SKIP:
        case cigarFlags::SOFT_CLIP:
        case cigarFlags::HARD_CLIP:
        case cigarFlags::PADDING:
        case cigarFlags::SEQ_MATCH:
        case cigarFlags::SEQ_MISMATCH:
            if ( length < cigarItemCacheSize and length >= 0 )
                return cigarItemCache[std::make_pair( flag, length )];
            else
                return rollCigarItem( flag, length );
        default:
            throw utils::wecall_exception( "Cigar flag number: " + std::to_string( static_cast< int >( flag ) ) +
                                            " not implemented." );
        }
    }

    std::shared_ptr< CigarItem > CigarItemBuilder::roll( const uint32_t & cigarBitFlag )
    {
        const auto flag = static_cast< cigarFlags >( cigarBitFlag & BAM_CIGAR_MASK );
        const auto length = cigarBitFlag >> BAM_CIGAR_SHIFT;
        return CigarItemBuilder::roll( flag, length );
    }

    std::shared_ptr< CigarItem > CigarItemBuilder::roll( const std::string & cigarItemString )
    {
        boost::regex regex( "^(\\d+)([MIDNSHP=X])$" );
        boost::cmatch what;
        if ( boost::regex_match( cigarItemString.c_str(), what, regex ) )
        {
            auto length = boost::lexical_cast< int64_t >( what[1] );
            char type = what[2].str()[0];
            cigarFlags flag;

            switch ( type )
            {
            case 'M':
                flag = cigarFlags::MATCH;
                break;
            case 'I':
                flag = cigarFlags::INSERTION;
                break;
            case 'D':
                flag = cigarFlags::DELETION;
                break;
            case 'N':
                flag = cigarFlags::SKIP;
                break;
            case 'S':
                flag = cigarFlags::SOFT_CLIP;
                break;
            case 'H':
                flag = cigarFlags::HARD_CLIP;
                break;
            case 'P':
                flag = cigarFlags::PADDING;
                break;
            case '=':
                flag = cigarFlags::SEQ_MATCH;
                break;
            case 'X':
                flag = cigarFlags::SEQ_MISMATCH;
                break;
            default:
                throw utils::wecall_exception( "Cigar flag type: " + std::string( 1, type ) + " not implemented." );
            }
            return roll( flag, length );
        }
        else
        {
            throw utils::wecall_exception( cigarItemString + " is not a valid cigar item" );
        }
    }

    Cigar::Cigar( const uint32_t nOperations, const uint32_t * cigarPtr )
    {
        m_cigarItems.resize( nOperations );
        for ( uint32_t i = 0; i < nOperations; ++i )
        {
            m_cigarItems[i] = CigarItemBuilder::roll( cigarPtr[i] );
        }
        m_lengthBeforeRefStartPos = this->computeLengthBeforeRefStartPos();
        m_lengthAfterRefEndPos = this->computeLengthAfterRefEndPos();
    }

    Cigar::Cigar( const std::string & cigarString )
    {
        boost::regex regex( "([\\d]+[MIDNSHP=X])" );
        auto cigar_begin = boost::sregex_iterator( cigarString.begin(), cigarString.end(), regex );
        auto cigar_end = boost::sregex_iterator();
        for ( boost::sregex_iterator i = cigar_begin; i != cigar_end; ++i )
        {
            auto cigarItem = CigarItemBuilder::roll( i->str() );
            m_cigarItems.push_back( cigarItem );
        }

        m_lengthBeforeRefStartPos = this->computeLengthBeforeRefStartPos();
        m_lengthAfterRefEndPos = this->computeLengthAfterRefEndPos();
    }

    int64_t Cigar::length() const
    {
        return std::accumulate( m_cigarItems.cbegin(), m_cigarItems.cend(), 0,
                                []( int64_t curSum, const std::shared_ptr< CigarItem > & item )
                                {
                                    return curSum + item->length();
                                } );
    }

    int64_t Cigar::lengthInRef() const
    {
        return std::accumulate( m_cigarItems.cbegin(), m_cigarItems.cend(), 0,
                                []( int64_t curSum, const std::shared_ptr< CigarItem > & item )
                                {
                                    return curSum + item->lengthInRef();
                                } );
    }

    int64_t Cigar::lengthInSeq() const
    {
        return std::accumulate( m_cigarItems.cbegin(), m_cigarItems.cend(), 0,
                                []( int64_t curSum, const std::shared_ptr< CigarItem > & item )
                                {
                                    return curSum + item->lengthInSeq();
                                } );
    }

    std::string Cigar::toString() const
    {
        return std::accumulate( m_cigarItems.cbegin(), m_cigarItems.cend(), std::string(),
                                []( std::string curMsg, const std::shared_ptr< CigarItem > & item )
                                {
                                    return curMsg + item->toString();
                                } );
    }

    int64_t Cigar::computeLengthBeforeRefStartPos() const
    {
        int64_t result( 0 );
        for ( auto it = m_cigarItems.cbegin(); it != m_cigarItems.cend() and ( *it )->lengthInRef() == 0; ++it )
        {
            result += ( *it )->length();
        }
        return result;
    }

    int64_t Cigar::computeLengthAfterRefEndPos() const
    {
        int64_t result( 0 );
        for ( auto it = m_cigarItems.crbegin(); it != m_cigarItems.crend() and ( *it )->lengthInRef() == 0; ++it )
        {
            result += ( *it )->length();
        }
        return result;
    }

    referencePositions_t Cigar::getRefPositions( int64_t startPos ) const
    {
        referencePositions_t refPositions;
        auto currentPos = startPos;
        for ( const auto item : m_cigarItems )
        {
            const auto itemsRefPositions = item->getRelativeRefPositions( currentPos );
            refPositions.insert( refPositions.end(), itemsRefPositions.begin(), itemsRefPositions.end() );
            currentPos += item->lengthInRef();
        }
        return refPositions;
    }

    int64_t Cigar::lengthInSeqWithoutSoftClipping() const
    {
        const auto accumulator = []( int64_t curSum, const std::shared_ptr< CigarItem > & item )
        {
            if ( item->isSoftClipped() )
            {
                return curSum;
            }
            else
            {
                return curSum + item->lengthInSeq();
            }
        };
        return std::accumulate( m_cigarItems.cbegin(), m_cigarItems.cend(), 0, accumulator );
    }

    std::pair< int64_t, int64_t > Cigar::stripSoftClipping()
    {
        int64_t lengthOfFront = 0L;
        int64_t lengthOfBack = 0L;

        if ( m_cigarItems.size() > 0 and m_cigarItems.front()->isSoftClipped() )
        {
            lengthOfFront = m_cigarItems.front()->lengthInSeq();
            m_cigarItems.erase( m_cigarItems.begin() );
        }

        if ( m_cigarItems.size() > 0 and m_cigarItems.back()->isSoftClipped() )
        {
            lengthOfBack = m_cigarItems.back()->lengthInSeq();
            m_cigarItems.pop_back();
        }

        return {lengthOfFront, lengthOfBack};
    }

    utils::Interval Cigar::getInverseInterval( const utils::Interval & input ) const
    {
        offsetsPtr_t currentPositions = std::make_shared< Offsets >( 0L, 0L );

        auto currentCigarItem = m_cigarItems.cbegin();
        for ( ; currentCigarItem != m_cigarItems.cend() and
                    currentPositions->ref + ( *currentCigarItem )->lengthInRef() < input.start();
              ++currentCigarItem )
        {
            ( *currentCigarItem )->moveOffsets( currentPositions );
        }

        int64_t inverseStart = currentPositions->read;
        if ( currentCigarItem != m_cigarItems.cend() and
             ( *currentCigarItem )->lengthInSeq() == ( *currentCigarItem )->lengthInRef() )
        {
            inverseStart += ( input.start() - currentPositions->ref );
        }

        for ( ; currentCigarItem != m_cigarItems.cend() and
                    currentPositions->ref + ( *currentCigarItem )->lengthInRef() <= input.end();
              ++currentCigarItem )
        {
            ( *currentCigarItem )->moveOffsets( currentPositions );
        }

        int64_t inverseEnd = currentPositions->read;
        if ( currentCigarItem != m_cigarItems.cend() and
             ( *currentCigarItem )->lengthInSeq() == ( *currentCigarItem )->lengthInRef() )
        {
            inverseEnd += ( input.end() - currentPositions->ref );
        }

        return utils::Interval( inverseStart, inverseEnd );
    }

}  // namespace alignment
}  // namespace wecall
