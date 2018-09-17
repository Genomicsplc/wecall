// All content Copyright (C) 2018 Genomics plc
#ifndef READ_CONTAINER_HPP
#define READ_CONTAINER_HPP

#include "io/readIntervalTree.hpp"
#include "utils/interval.hpp"

#include <utility>
#include <map>
#include <string>
#include <cstdint>

namespace echidna
{
namespace io
{
    using readRange_t = std::pair< readIt_t, readIt_t >;

    /// Return a pair of iterators spanning the specified interval
    ///
    class RegionsReads
    {
        // Class to wrap iteration of reads in weCall.
        class iterator : public std::iterator< std::forward_iterator_tag, io::Read >
        {
        public:
            iterator( readIt_t current, const RegionsReads * parent );

            const io::Read & operator*() const { return m_current.operator*(); }
            io::Read & operator*() { return m_current.operator*(); }
            io::readPtr_t getSharedPtr() { return m_current.operator->(); }

            iterator & operator++();
            bool operator!=( iterator rhs ) const { return this->m_current != rhs.m_current; }
            bool operator==( iterator rhs ) const { return this->m_current == rhs.m_current; }

        private:
            void moveToValidInterval();
            readIt_t m_current;
            const RegionsReads * m_parent;
        };

    public:
        RegionsReads( const caller::SetRegions & regions, const readRange_t & reads, phred_t minMappingQuality );

        iterator begin() const;
        iterator end() const;

        RegionsReads getSubRegionReads( const caller::SetRegions & subRegions ) const;

        caller::SetRegions getRegions() const { return m_regions; }

        std::string toString() const;

    private:
        bool isValid( const iterator & it ) const
        {
            if ( it == this->end() )
            {
                return true;
            }
            else if ( m_regions.empty() )
            {
                return false;
            }
            else if ( it.operator*().getMappingQuality() < m_minMappingQuality )
            {
                return false;
            }

            const utils::Interval interval = it.operator*().getMaximalReadInterval();
            return m_regions.overlaps( caller::Region( m_regions.cbegin()->contig(), interval ) );
        }

        const caller::SetRegions m_regions;
        const readRange_t m_reads;
        const phred_t m_minMappingQuality;
    };

    typedef std::map< std::string, RegionsReads > perSampleRegionsReads_t;

    perSampleRegionsReads_t reduceRegionSet( const perSampleRegionsReads_t & reads,
                                             const caller::SetRegions & regions );

    int64_t minReadStartPos( const RegionsReads & readRange );
    int64_t maxReadAlignedPos( const RegionsReads & readRange );
    std::pair< int64_t, int64_t > readsAlignedStartEnd( const RegionsReads & readRange );

    int64_t perSampleMaxAlignedReadLength( const perSampleRegionsReads_t & perSampleReadRanges );

    int64_t maxAlignedReadLength( const RegionsReads & readRange );
    int64_t maxReadLength( const RegionsReads & readRange );
    int64_t perSampleMaxReadLength( const perSampleRegionsReads_t & perSampleReadRanges );

    int64_t maxReadCigarLength( const RegionsReads & readRange );
    int64_t perSampleMaxReadCigarLength( const perSampleRegionsReads_t & perSampleReadRanges );

    bool readsOverlappingRegions( const perSampleRegionsReads_t & perSampleReads,
                                  const caller::Region & region1,
                                  const caller::Region & region2 );
}
}

#endif
