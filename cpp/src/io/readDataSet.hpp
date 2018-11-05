// All content Copyright (C) 2018 Genomics plc
#ifndef READ_DATASET_HPP
#define READ_DATASET_HPP

#include "common.hpp"
#include "io/readRange.hpp"

#include "io/readfilters/readFilterAndTrimmer.hpp"
#include "caller/params.hpp"
#include "caller/region.hpp"
#include "utils/logging.hpp"
#include "utils/interval.hpp"

namespace wecall
{
namespace io
{
    using readData_t = std::map< std::string, readIntervalTree_t >;
    /// Stores reads from >= 1 samples.
    class ReadDataset
    {
    public:
        ReadDataset( std::vector< std::string > sampleNames, caller::Region region );

        /// Disabled copy constructor. No copying of datasets
        ReadDataset( const ReadDataset & rhs ) = delete;

        /// Disabled assignment operator. No copying of datasets
        ReadDataset & operator=( const ReadDataset & ) = delete;

        /// Returns the set of all samples observed in the loaded data
        std::vector< std::string > getSampleNames() const { return m_samples; }

        caller::Region region() const { return m_region; }

        void updateRegionEnd( int64_t newRegionEnd )
        {
            m_region = caller::Region( m_region.contig(), m_region.start(), newRegionEnd );
        }

        perSampleRegionsReads_t getRegionsReads( const caller::SetRegions & setRegions,
                                                 phred_t minMappingQuality ) const;

        /// Returns a map of between sample names and read ranges.
        ///
        perSampleRegionsReads_t getAllReads( phred_t minMappingQuality ) const;

        bool isEmpty() const { return m_empty; }

        /// Inserts per sample the new Read into the data set. Except if the read is
        // flagged unmapped and/or getStartPos == getAlignedPos.
        void insertRead( const std::string & sampleName, readPtr_t readPtr );

    private:
        caller::Region m_region;
        std::vector< std::string > m_samples;
        readData_t m_intervalTreeData;
        bool m_empty;
    };

    using readDataset_t = std::shared_ptr< const ReadDataset >;
}
}

#endif
