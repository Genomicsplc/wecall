// All content Copyright (C) 2018 Genomics plc
#include <vector>
#include <set>

#include <boost/algorithm/string/join.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include "utils/partition.hpp"
#include "utils/identity.hpp"
#include "caller/region.hpp"
#include "vcf/reader.hpp"
#include "caller/regionUtils.hpp"
#include "io/fastaFile.hpp"
#include "io/bedFile.hpp"

namespace echidna
{
namespace caller
{
    partitionedRegions_t DataRegionsBuilderBase::build()
    {
        const auto inputRegions = this->getRegions();

        auto regions = this->getValidRegions( inputRegions );

        if ( not std::is_sorted( regions.begin(), regions.end(), RegionComp( m_fastaIndexFile ) ) )
        {
            ECHIDNA_LOG( WARNING, "Regions specified are not in order. Sorting them by contig:start-end." );

            std::sort( regions.begin(), regions.end(), RegionComp( m_fastaIndexFile ) );
            ECHIDNA_LOG( SUPER_DEBUG,
                         "Regions sorted to: " + boost::algorithm::join( getRegionStrings( regions ), "," ) );
        }

        auto sameContig = []( const Region & a, const Region & b ) -> bool
        {
            return a.contig() == b.contig();
        };
        auto partitionedRegionsUnMerged = utils::functional::partition( regions, sameContig );

        partitionedRegions_t partitionedRegions;

        for ( const auto & partitionUnMerged : partitionedRegionsUnMerged )
        {
            partitionedRegions.emplace_back( mergeRegions( partitionUnMerged ) );
        }

        return partitionedRegions;
    }

    regions_t DataRegionsBuilderBase::getValidRegions( const regions_t & regions ) const
    {
        // Check the regions
        const auto allowedContigs = m_fastaIndexFile.contigs();

        regions_t goodRegions;
        std::set< std::string > badContigs;

        for ( const auto & region : regions )
        {
            const auto it = allowedContigs.find( region.contig() );

            if ( it == allowedContigs.cend() )
            {
                badContigs.insert( region.contig() );
            }
            else
            {
                goodRegions.emplace_back( region );
            }
        }

        if ( not badContigs.empty() )
        {
            ECHIDNA_LOG( WARNING, "Contig(s) " + boost::algorithm::join( badContigs, "," ) +
                                      " are not contained in reference." );
        }

        return goodRegions;
    }

    regions_t DataRegionsBuilder::getRegions() const
    {
        regions_t regions;
        if ( m_inputRegionStrings.empty() )
        {
            regions = this->getDefaultRegions();
        }
        else
        {
            std::vector< bool > areStringsBedFiles;
            std::transform( m_inputRegionStrings.begin(), m_inputRegionStrings.end(),
                            std::back_inserter( areStringsBedFiles ), this->isBedFile );

            using echidna::utils::functional::identity;

            if ( std::all_of(areStringsBedFiles.begin(), areStringsBedFiles.end(), identity< bool >))
            {
                regions = this->getRegionsFromBedFiles();
            }
            else if ( std::none_of(areStringsBedFiles.begin(), areStringsBedFiles.end(), identity< bool >))
            {
                regions = this->getRegionsFromStrings();
            }
            else
            {
                throw utils::echidna_exception( "Can not have mixture of BED files and region strings" );
            }
        }
        return regions;
    }

    regions_t DataRegionsBuilder::getDefaultRegions() const
    {
        regions_t regions;
        const auto contigs = m_fastaIndexFile.contigs();

        for ( const auto & contigName : m_fastaIndexFile.standardContigs() )
        {
            regions.emplace_back( contigName, contigs.at( contigName ) );
        }

        ECHIDNA_ASSERT( std::is_sorted( regions.begin(), regions.end(), RegionComp( m_fastaIndexFile ) ),
                        "standard contigs in FastaFile not sorted" );

        return regions;
    }

    regions_t DataRegionsBuilder::getRegionsFromBedFiles() const
    {
        regions_t regions;
        for ( const auto & inputRegionString : m_inputRegionStrings )
        {
            auto newRegions = io::BedFile( inputRegionString ).getRegions();
            regions.insert( regions.end(), newRegions.begin(), newRegions.end() );
        }
        return regions;
    }

    regions_t DataRegionsBuilder::getRegionsFromStrings() const
    {
        regions_t regions;
        const auto contigs = m_fastaIndexFile.contigs();
        for ( auto regionString : m_inputRegionStrings )
        {
            regions.push_back( parseRegionString( regionString, contigs ) );
        }

        return regions;
    }

    bool DataRegionsBuilder::isBedFile( std::string inputString )
    {
        auto possPath = boost::filesystem::path( inputString );
        bool isBed = possPath.extension().string() == ".bed";
        bool isZippedBed = possPath.extension().string() == ".gz" and possPath.stem().extension().string() == ".bed";
        return isBed or isZippedBed;
    }

    std::vector< std::string > getRegionStrings( const std::vector< Region > & regions )
    {
        std::vector< std::string > regionStrings;
        for ( const auto & region : regions )
        {
            regionStrings.emplace_back( region.toString() );
        }
        return regionStrings;
    }

    std::vector< Region > getContigsFromRegions(
        const std::vector< Region > & regions,
        const std::map< std::string, utils::Interval > & contigsFromFastaIndexFile )
    {
        std::vector< Region > contigs;
        std::set< std::string > contigStrings;
        for ( const auto & region : regions )
        {
            if ( contigStrings.find( region.contig() ) == contigStrings.end() )
            {
                contigStrings.insert( region.contig() );
                contigs.emplace_back( region.contig(), contigsFromFastaIndexFile.at( region.contig() ) );
            }
        }
        return contigs;
    }

    boost::regex chromOnly( "^" + chromMatch + "$" );
    boost::regex fullRegion( "^" + chromMatch + ":" + posMatch + "-" + posMatch + "$" );
    Region parseRegionString( std::string regionString, const std::map< std::string, utils::Interval > & contigs )
    {
        boost::cmatch what;
        if ( boost::regex_match( regionString.c_str(), what, chromOnly ) )
        {
            std::string contig = what[1].str();

            // No interval provided, so find it in the reference
            auto it = contigs.find( contig );
            if ( it != contigs.end() )
            {
                return Region( what[1].str(), it->second );
            }
            else
            {
                return Region( what[1].str(), 0, 0 );
            }
        }
        else if ( boost::regex_match( regionString.c_str(), what, fullRegion ) )
        {
            return Region( what[1].str(), utils::Interval( boost::lexical_cast< int64_t >( what[2] ),
                                                           boost::lexical_cast< int64_t >( what[3] ) ) );
        }
        else
        {
            throw utils::echidna_exception( "Invalid region specification format: " + regionString +
                                            ". Expecting: chrom or chrom:start-end" );
        }
    }

    std::vector< Region > mergeRegions( const std::vector< Region > & regions )
    {
        if ( regions.empty() )
        {
            return {};
        }

        const auto contig = regions.front().contig();
        const auto sameContig = [contig]( const Region & region )
        {
            return region.contig() == contig;
        };

        ECHIDNA_ASSERT( std::all_of( regions.begin(), regions.end(), sameContig ),
                        "Can not merge regions on different contigs" );

        ECHIDNA_ASSERT( std::is_sorted( regions.begin(), regions.end(), RegionIntervalComp() ),
                        "Merge Regions requires sorted regions" );

        std::vector< Region > output = {regions.front()};
        for ( const auto & region : regions )
        {
            if ( output.back().overlaps( region ) )
            {
                output.back().combine( region );
            }
            else if ( output.back().end() == region.start() )
            {
                output.back().combine( region );
            }
            else
            {
                output.push_back( region );
            }
        }

        return output;
    }

    regions_t padRegions( const regions_t & regions, const int64_t paddingAmount )
    {
        if ( regions.empty() )
        {
            return {};
        }
        regions_t paddedRegions;
        paddedRegions.reserve( regions.size() );
        const auto pad = [paddingAmount]( const Region & unpadded )
        {
            return unpadded.getPadded( paddingAmount );
        };

        std::transform( regions.begin(), regions.end(), std::back_inserter( paddedRegions ), pad );
        return mergeRegions( paddedRegions );
    }

    partitionedRegions_t padPartitionedRegions( const partitionedRegions_t & partitionedRegions,
                                                const int64_t paddingAmount )
    {
        if ( partitionedRegions.empty() )
        {
            return {};
        }
        partitionedRegions_t paddedPartitionedRegions;
        paddedPartitionedRegions.reserve( partitionedRegions.size() );

        const auto pad = [paddingAmount]( const regions_t & input )
        {
            return padRegions( input, paddingAmount );
        };

        std::transform( partitionedRegions.begin(), partitionedRegions.end(),
                        std::back_inserter( paddedPartitionedRegions ), pad );
        return paddedPartitionedRegions;
    }

    bool RegionComp::operator()( const Region & a, const Region & b )
    {
        if ( m_fastaIndexFile.contigStart( a.contig() ) == m_fastaIndexFile.contigStart( b.contig() ) )
        {
            return RegionIntervalComp()( a, b );
        }

        return m_fastaIndexFile.contigStart( a.contig() ) < m_fastaIndexFile.contigStart( b.contig() );
    }
}
}
