// All content Copyright (C) 2018 Genomics plc
#include "io/readRange.hpp"
#include "alignment/cigarItems.hpp"

namespace wecall
{
namespace io
{

    perSampleRegionsReads_t reduceRegionSet( const perSampleRegionsReads_t & reads, const caller::SetRegions & regions )
    {
        perSampleRegionsReads_t ret;
        for ( const auto & pair : reads )
        {
            ret.emplace( std::piecewise_construct, std::forward_as_tuple( pair.first ),
                         std::forward_as_tuple( pair.second.getSubRegionReads( regions ) ) );
        }
        return ret;
    }

    RegionsReads RegionsReads::getSubRegionReads( const caller::SetRegions & subRegions ) const
    {
        const auto subRegionSpan = subRegions.getSpan();
        WECALL_ASSERT( m_regions.getSpan().contains( subRegionSpan ),
                        "Span of sub regions ( " + subRegions.toString() +
                            " ) required to be contained in m_regions (" + m_regions.toString() + ")" );

        readRange_t readRange = m_reads;
        readRange.first.moveToBeginOfSubrange( subRegionSpan.interval() );
        return RegionsReads( subRegions, readRange, m_minMappingQuality );
    }

    //-----------------------------------------------------------------------------------------

    int64_t minReadStartPos( const RegionsReads & readRange )
    {
        if ( readRange.end() == readRange.begin() )
        {
            return alignment::noPos;
        }
        else
        {
            auto comp = []( const io::Read & left, const io::Read & right )
            {
                return left.getStartPos() < right.getStartPos();
            };
            auto min_element = std::min_element( readRange.begin(), readRange.end(), comp );
            return ( *min_element ).getStartPos();
        }
    }

    //-----------------------------------------------------------------------------------------

    int64_t maxReadAlignedPos( const RegionsReads & readRange )
    {
        if ( readRange.end() == readRange.begin() )
        {
            return alignment::noPos;
        }
        else
        {
            auto comp = []( const io::Read & left, const io::Read & right )
            {
                return left.getAlignedEndPos() < right.getAlignedEndPos();
            };
            // Note odd implementation of ReadEndPosComp which returns > not <.
            auto max_element = std::max_element( readRange.begin(), readRange.end(), comp );
            return ( *max_element ).getAlignedEndPos();
        }
    }

    //-----------------------------------------------------------------------------------------

    std::pair< int64_t, int64_t > readsAlignedStartEnd( const RegionsReads & readRange )
    {

        const auto lastReadEnd = maxReadAlignedPos( readRange );
        const auto firstReadStart = minReadStartPos( readRange );

        return std::make_pair( firstReadStart, lastReadEnd );
    }

    //-----------------------------------------------------------------------------------------

    int64_t maxAlignedReadLength( const RegionsReads & readRange )
    {
        int64_t maxReadLength = 0;  // so the std::max compiles.
        for ( const auto & read : readRange )
        {
            maxReadLength = std::max( maxReadLength, read.getAlignedLength() );
        }
        return maxReadLength;
    }

    //-----------------------------------------------------------------------------------------

    int64_t perSampleMaxAlignedReadLength( const perSampleRegionsReads_t & perSampleReadRanges )
    {
        int64_t maxReadLength = 0;

        for ( auto sampleReadsIt = perSampleReadRanges.begin(); sampleReadsIt != perSampleReadRanges.end();
              ++sampleReadsIt )
        {
            const auto thisItMax = maxAlignedReadLength( sampleReadsIt->second );
            maxReadLength = std::max( maxReadLength, thisItMax );
        }
        return maxReadLength;
    }

    int64_t maxReadLength( const RegionsReads & readRange )
    {
        int64_t maxLength = 0;  // so the std::max compiles.
        for ( const auto & read : readRange )
        {
            maxLength = std::max( maxLength, read.getLength() );
        }
        return maxLength;
    }

    int64_t maxReadCigarLength( const RegionsReads & readRange )
    {
        int64_t maxLength = 0;  // so the std::max compiles.
        for ( const auto & read : readRange )
        {
            maxLength = std::max( maxLength, read.cigar().length() );
        }
        return maxLength;
    }

    int64_t perSampleMaxReadLength( const perSampleRegionsReads_t & perSampleReadRanges )
    {
        int64_t maxLength = 0;

        for ( auto sampleReadsIt = perSampleReadRanges.begin(); sampleReadsIt != perSampleReadRanges.end();
              ++sampleReadsIt )
        {
            const auto thisItMax = maxReadLength( sampleReadsIt->second );
            maxLength = std::max( maxLength, thisItMax );
        }
        return maxLength;
    }

    int64_t perSampleMaxReadCigarLength( const perSampleRegionsReads_t & perSampleReadRanges )
    {
        int64_t maxLength = 0;

        for ( auto sampleReadsIt = perSampleReadRanges.begin(); sampleReadsIt != perSampleReadRanges.end();
              ++sampleReadsIt )
        {
            const auto thisItMax = maxReadCigarLength( sampleReadsIt->second );
            maxLength = std::max( maxLength, thisItMax );
        }
        return maxLength;
    }

    RegionsReads::RegionsReads( const caller::SetRegions & regions,
                                const readRange_t & reads,
                                phred_t minMappingQuality )
        : m_regions( regions ), m_reads( reads ), m_minMappingQuality( minMappingQuality )
    {
        WECALL_ASSERT( m_regions.allSameContig(), "Current only deal with regions of one contig at a time." );
    }

    RegionsReads::iterator::iterator( readIt_t current, const RegionsReads * parent )
        : m_current( current ), m_parent( parent )
    {
    }

    void RegionsReads::iterator::moveToValidInterval()
    {
        while ( not m_parent->isValid( *this ) )
        {
            ++m_current;
        }
    }

    std::string RegionsReads::toString() const
    {
        std::stringstream sstr;
        sstr << "RegionReads for:\t" << m_regions.toString();
        return sstr.str();
    }

    RegionsReads::iterator RegionsReads::begin() const
    {
        iterator iteratorBegin( m_reads.first, this );
        if ( not this->isValid( iteratorBegin ) )
        {
            ++iteratorBegin;
        }
        return iteratorBegin;
    }

    RegionsReads::iterator RegionsReads::end() const { return iterator( m_reads.second, this ); }

    RegionsReads::iterator & RegionsReads::iterator::operator++()
    {
        ++m_current;
        this->moveToValidInterval();
        return *this;
    }

    //-----------------------------------------------------------------------------------------

    bool readsOverlappingRegions( const perSampleRegionsReads_t & perSampleReads,
                                  const caller::Region & region1,
                                  const caller::Region & region2 )
    {
        // function includes read pairs as they share the same QName
        perSampleRegionsReads_t ret;
        std::vector< std::pair< std::string, int64_t > > read2Ids;
        for ( const auto & pair : perSampleReads )
        {
            const io::RegionsReads reads1 = pair.second.getSubRegionReads( region1 );
            const io::RegionsReads reads2 = pair.second.getSubRegionReads( region2 );

            read2Ids.clear();
            for ( const auto & read : reads2 )
            {
                // the qName is shared between read pair. This is to make it a unique id
                read2Ids.push_back( std::make_pair( read.getQName(), read.getInsertSize() ) );
            }
            std::sort( read2Ids.begin(), read2Ids.end() );

            for ( const auto & read : reads1 )
            {
                auto found = std::binary_search( read2Ids.cbegin(), read2Ids.cend(),
                                                 std::make_pair( read.getQName(), read.getInsertSize() ) );
                if ( found )
                    return true;
            }
        }
        return false;
    }
}
}
