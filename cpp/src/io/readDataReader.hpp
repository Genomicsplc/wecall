// All content Copyright (C) 2018 Genomics plc
#ifndef READ_DATAREADER_HPP
#define READ_DATAREADER_HPP

#include "common.hpp"

#include "io/bamFile.hpp"
#include "io/readDataSet.hpp"
#include "caller/params.hpp"
#include "utils/logging.hpp"
#include "caller/region.hpp"

namespace echidna
{
namespace io
{
    using readMap_t = std::map< std::string, std::vector< readPtr_t > >;
    /// Loads reads from multiple sources.
    class ReadDataReader
    {
    public:
        class BlockIterator
        {
        public:
            BlockIterator( const ReadDataReader * reader,
                           const caller::Region & region,
                           utils::referenceSequencePtr_t refSequence )
                : m_reader( reader ),
                  m_refSequence( refSequence ),
                  m_region( region ),
                  m_curPos( region.start() ),
                  m_memUsed( 0 )
            {
            }

            /// Read in the next block of reads
            ///
            /// @return read data set or NULL if no more reads.
            readDataset_t getReadDatasetForNextBlock();
            bool isLastBlock();
            void chopCurrentBlock( int64_t prematureBlockEnd );

        private:
            readMap_t takeBite( std::vector< bamFileIteratorPtr_t > iterators, int64_t biteToPos );
            bool isFull() { return m_memUsed > m_reader->m_memLimit; }
            bool isAlmostFull() { return m_memUsed > ( m_reader->m_memLimit * 0.8 ) and not this->isFull(); }

        private:
            const ReadDataReader * m_reader;
            utils::referenceSequencePtr_t m_refSequence;
            const caller::Region m_region;
            int64_t m_curPos;
            int64_t m_memUsed;
        };

        ReadDataReader( const caller::params::System & systemParams,
                        const caller::params::PrivateSystem & privateSystemParams,
                        const caller::params::Filters & filterParams,
                        const int64_t biteSize,
                        const std::vector< std::string > & bamFiles )
            : m_maxBlockSize( systemParams.m_maxBlockSize ),
              m_biteSize( biteSize ),
              m_memLimit( privateSystemParams.m_memLimit * 1024 * 1024 ),
              m_readFilterAndTrimmer( filterParams )
        {
            initDataSources( bamFiles );
        }

        ReadDataReader( int64_t maxBlockSize,
                        int64_t memLimitBytes,
                        const caller::params::Filters & filterParams,
                        const int64_t biteSize,
                        const std::vector< std::string > & dataSources )
            : m_maxBlockSize( maxBlockSize ),
              m_biteSize( biteSize ),
              m_memLimit( memLimitBytes ),
              m_readFilterAndTrimmer( filterParams )
        {
            initDataSources( dataSources );
        }

        /// Disable copying as instances will be holding open data sources
        ReadDataReader( const ReadDataReader & rhs ) = delete;

        /// Disable assignment as instances will be holding open data sources
        ReadDataReader & operator=( const ReadDataReader & ) = delete;

        /// Returns an iterator to read from the specified region
        ///
        /// @param region The region to read
        ReadDataReader::BlockIterator readRegion( const caller::Region & region,
                                                  utils::referenceSequencePtr_t refSequence );

        /// Returns the set of all samples observed in the loaded data
        std::vector< std::string > getSampleNames() const { return m_samples; }

        void initDataSources( const std::vector< std::string > & dataSources );

    private:
        /// Add a data source to this reader. All datasources that have been added will be used to load data.
        void addDataSource( bamFilePtr_t dataSource );

    private:
        std::vector< bamFilePtr_t > m_dataSources;
        std::vector< std::string > m_samples;
        int64_t m_maxBlockSize;
        int64_t m_biteSize;
        int64_t m_memLimit;

        ReadFilterAndTrimmer m_readFilterAndTrimmer;
    };
}
}

#endif
