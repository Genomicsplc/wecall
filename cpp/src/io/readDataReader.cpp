// All content Copyright (C) 2018 Genomics plc
#include <chrono>
#include "utils/timer.hpp"
#include "io/readDataReader.hpp"

namespace wecall
{
namespace io
{

    //-----------------------------------------------------------------------------------------

    void ReadDataReader::initDataSources( const std::vector< std::string > & dataSources )
    {
        for ( const std::string & sourceName : dataSources )
        {
            WECALL_LOG( DEBUG, "Adding data source " << sourceName << " to read dataset" );
            addDataSource( std::make_shared< BamFile >( sourceName ) );
        }
    }

    //-----------------------------------------------------------------------------------------

    void ReadDataReader::addDataSource( bamFilePtr_t dataSource )
    {
        m_dataSources.push_back( dataSource );

        for ( const auto & sampleName : dataSource->getSampleNames() )
        {
            if ( std::find( m_samples.begin(), m_samples.end(), sampleName ) == m_samples.end() )
            {
                WECALL_LOG( DEBUG, "Found sample " << sampleName << " in input data" );
                m_samples.push_back( sampleName );
            }

            else
            {
                WECALL_LOG( WARNING, "Found duplicate sample " << sampleName << " in input data" );
            }
        }
    }

    //-----------------------------------------------------------------------------------------

    ReadDataReader::BlockIterator ReadDataReader::readRegion( const caller::Region & region,
                                                              utils::referenceSequencePtr_t refSequence )
    {
        return BlockIterator( this, region, refSequence );
    }

    //-----------------------------------------------------------------------------------------

    readDataset_t ReadDataReader::BlockIterator::getReadDatasetForNextBlock()
    {
        if ( m_curPos >= m_region.end() )
        {
            return nullptr;
        }

        m_memUsed = 0;

        int64_t blockStart = m_curPos;
        int64_t nBlocksRemaining = ( ( m_region.end() - m_curPos - 1 ) / m_reader->m_maxBlockSize ) + 1;
        int64_t blockEnd = m_curPos + ( m_region.end() - m_curPos ) / nBlocksRemaining;

        const caller::Region blockRegion( m_region.contig(), blockStart, blockEnd );
        auto dataset = std::make_shared< ReadDataset >( m_reader->getSampleNames(), blockRegion );

        std::vector< bamFileIteratorPtr_t > iterators;
        const auto paddedBlockRegion = blockRegion.getPadded( constants::bamFetchRegionPadding );
        for ( auto dataSource : m_reader->m_dataSources )
        {
            auto bamFileIterator = dataSource->readRegion( paddedBlockRegion, m_refSequence );
            if ( bamFileIterator != nullptr )
            {
                iterators.push_back( bamFileIterator );
            }
        }

        while ( m_curPos < blockEnd )
        {
            auto biteToPos = std::min( blockEnd, m_curPos + m_reader->m_biteSize );
            readMap_t readData = this->takeBite( iterators, biteToPos );

            // If we reached full taking this bite, discard it and consider the block complete after last bite.
            if ( this->isFull() )
            {
                break;
            }

            for ( const auto & sampleReads : readData )
            {
                //                    filtered_reads = m_readFilterAndTrimmer.filter(sampleReads.second);
                for ( const auto read : sampleReads.second )
                {
                    if ( m_reader->m_readFilterAndTrimmer.trimAndFilter( read ) )
                    {
                        dataset->insertRead( sampleReads.first, read );
                    }
                }
            }
            m_curPos = biteToPos;

            // Special case - continue if some capacity and only a wafer thin mint left to digest!
            if ( this->isAlmostFull() and ( m_region.end() - m_curPos ) > m_reader->m_biteSize )
            {
                break;
            }
        }

        if ( m_curPos == blockStart )
        {
            // Couldn't manage a single bite - skip block and log a WARNING.
            m_curPos = std::min( blockEnd, m_curPos + m_reader->m_biteSize );
            blockEnd = m_curPos;
            WECALL_LOG( WARNING, "Skipping region " << caller::Region( m_region.contig(), blockStart, m_curPos )
                                                     << " due to exceptionally high coverage" );
        }
        else if ( m_curPos < blockEnd )
        {
            // Couldn't finish the whole meal - redefine block and log this INFO.
            blockEnd = m_curPos;
            WECALL_LOG( DEBUG, "Reducing block size due to high coverage in this region" );
        }

        dataset->updateRegionEnd( blockEnd );
        return dataset;
    }

    //-----------------------------------------------------------------------------------------

    readMap_t ReadDataReader::BlockIterator::takeBite( std::vector< bamFileIteratorPtr_t > iterators,
                                                       int64_t biteToPos )
    {
        readMap_t readMap;
        for ( auto iterator : iterators )
        {
            auto scopedTimerTrigger = iterator->timer();
            for ( ; iterator->hasReadData(); iterator->next() )
            {
                auto sampleRead = iterator->getReadData();

                std::string sampleName = sampleRead.first;
                readPtr_t read = sampleRead.second;
                if ( read->getStartPos() > biteToPos )
                {
                    // Gone past the end of the current bite for this sample move on to the next
                    break;
                }

                // Rough approximation
                m_memUsed += ( ( read->getLength() * 2 ) + 100 );
                if ( this->isFull() )
                {
                    return readMap;
                }

                readMap[sampleName].push_back( read );
            }
        }

        return readMap;
    }

    //-----------------------------------------------------------------------------------------

    bool ReadDataReader::BlockIterator::isLastBlock() { return ( m_curPos >= m_region.end() ); }

    //-----------------------------------------------------------------------------------------

    void ReadDataReader::BlockIterator::chopCurrentBlock( int64_t prematureBlockEnd )
    {
        WECALL_LOG( DEBUG, "Chopping current block at position: " << prematureBlockEnd
                                                                   << " to avoid breaking up a cluster of variants" );
        m_curPos = prematureBlockEnd;
    }
}
}
