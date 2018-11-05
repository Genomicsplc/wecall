// All content Copyright (C) 2018 Genomics plc
#ifndef CALLER_REGIONUTILS_HPP
#define CALLER_REGIONUTILS_HPP

#include "caller/region.hpp"
#include "io/fastaFile.hpp"

namespace wecall
{
namespace caller
{
    using partitionedRegions_t = std::vector< regions_t >;

    Region parseRegionString( std::string regionString, const std::map< std::string, utils::Interval > & contigs );

    class DataRegionsBuilderBase
    {
    public:
        DataRegionsBuilderBase( const io::FastaIndex & fastaIndex ) : m_fastaIndexFile( fastaIndex ) {}

        partitionedRegions_t build();

    protected:
        virtual regions_t getRegions() const = 0;
        regions_t getValidRegions( const regions_t & regions ) const;
        const io::FastaIndex & m_fastaIndexFile;
    };

    class DataRegionsBuilder : public DataRegionsBuilderBase
    {
    public:
        DataRegionsBuilder( std::vector< std::string > inputRegionStrings, const io::FastaIndex & fastaIndex )
            : DataRegionsBuilderBase( fastaIndex ), m_inputRegionStrings( inputRegionStrings )
        {
        }

    private:
        regions_t getDefaultRegions() const;
        regions_t getRegionsFromBedFiles() const;
        regions_t getRegionsFromStrings() const;

        static bool isBedFile( std::string inputString );

        regions_t getRegions() const override;

    private:
        std::vector< std::string > m_inputRegionStrings;
    };

    std::vector< std::string > getRegionStrings( const std::vector< Region > & regions );

    std::vector< Region > getContigsFromRegions(
        const std::vector< Region > & regions,
        const std::map< std::string, utils::Interval > & contigsFromFastaIndexFile );

    // merge sorted regions on the same contig
    regions_t mergeRegions( const regions_t & regions );

    regions_t padRegions( const regions_t & regions, const int64_t paddingAmount );
    partitionedRegions_t padPartitionedRegions( const partitionedRegions_t & partitionedRegions,
                                                const int64_t paddingAmount );

    struct RegionIntervalComp
    {
        bool operator()( const Region & a, const Region & b ) { return a.interval() < b.interval(); }
    };

    struct RegionComp
    {
        RegionComp( const io::FastaIndex & fastaIndexFile ) : m_fastaIndexFile( fastaIndexFile ) {}

        bool operator()( const Region & a, const Region & b );

        const io::FastaIndex & m_fastaIndexFile;
    };
}
}

#endif
