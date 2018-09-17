// All content Copyright (C) 2018 Genomics plc
#include "io/readDataSet.hpp"

namespace echidna
{
namespace io
{
    ReadDataset::ReadDataset( std::vector< std::string > sampleNames, caller::Region region )
        : m_region( region ), m_samples( sampleNames ), m_intervalTreeData(), m_empty( true )
    {
        for ( auto sampleName : m_samples )
        {
            m_intervalTreeData.emplace( std::piecewise_construct, std::forward_as_tuple( sampleName ),
                                        std::forward_as_tuple( region.start(), region.end() ) );
        }
    }

    perSampleRegionsReads_t ReadDataset::getRegionsReads( const caller::SetRegions & setRegions,
                                                          phred_t minMappingQuality ) const
    {
        const auto span = setRegions.getSpan();

        io::perSampleRegionsReads_t regionReads;
        for ( const auto & sampleTreePair : m_intervalTreeData )
        {
            regionReads.emplace(
                std::piecewise_construct, std::forward_as_tuple( sampleTreePair.first ),
                std::forward_as_tuple( setRegions, sampleTreePair.second.getSubRange( span.interval() ),
                                       minMappingQuality ) );
        }
        return regionReads;
    }

    //-----------------------------------------------------------------------------------------

    perSampleRegionsReads_t ReadDataset::getAllReads( phred_t minMappingQuality ) const
    {
        perSampleRegionsReads_t readRanges;
        caller::SetRegions setRegions( m_region );

        for ( auto & sampleTreePair : m_intervalTreeData )
        {
            readRanges.emplace(
                std::piecewise_construct, std::forward_as_tuple( sampleTreePair.first ),
                std::forward_as_tuple( setRegions, sampleTreePair.second.getFullRange(), minMappingQuality ) );
        }

        return readRanges;
    }

    //-----------------------------------------------------------------------------------------

    void ReadDataset::insertRead( const std::string & sampleName, readPtr_t readPtr )
    {
        m_empty = false;
        m_intervalTreeData.at( sampleName ).insert( readPtr );
    }

    //-----------------------------------------------------------------------------------------
}
}
