// All content Copyright (C) 2018 Genomics plc
#include <algorithm>
#include "variant/variantContainer.hpp"
#include "io/readDataSet.hpp"
#include "variant/variantGenerator.hpp"

namespace wecall
{
namespace variant
{

    //-----------------------------------------------------------------------------------------
    // VariantContainer::VariantCounts
    //-----------------------------------------------------------------------------------------

    //-----------------------------------------------------------------------------------------

    int64_t VariantContainer::VariantCounts::getPercentVariantCoverage() const
    {
        if ( m_totalReads == 0 )
        {
            return 0;
        }

        else
        {
            const double num = static_cast< double >( m_totalVariantSupportingReads );
            const auto denom = static_cast< double >( m_totalReads );
            const auto answer = ( 100.0 * num ) / denom;
            return static_cast< int64_t >( std::round( answer ) );
        }
    }

    //-----------------------------------------------------------------------------------------
    // VariantContainer
    //-----------------------------------------------------------------------------------------

    //-----------------------------------------------------------------------------------------

    void VariantContainer::addCandidateVariant( const varPtr_t & candidate, const double prior )
    {
        const auto findResult = m_variantCountsPerSample.find( candidate );
        varPtr_t var;
        if ( findResult == m_variantCountsPerSample.end() )
        {
            m_variantCountsPerSample[candidate] = std::map< std::string, VariantCounts >();
            var = candidate;
        }
        else
        {
            var = findResult->first;
        }

        var->disableFiltering();
        var->prior( prior );
    }

    void VariantContainer::addGenotypingVariant( const varPtr_t & variant )
    {
        const auto findResult = m_variantCountsPerSample.find( variant );
        varPtr_t var;
        if ( findResult == m_variantCountsPerSample.end() )
        {
            m_variantCountsPerSample[variant] = std::map< std::string, VariantCounts >();
            var = variant;
        }
        else
        {
            var = findResult->first;
        }

        var->disableFiltering();
        var->setGenotypingVariant();
    }

    void VariantContainer::addVariantsFromRead( io::readPtr_t readPtr,
                                                const std::vector< varPtr_t > & variants,
                                                const std::vector< breakpointPtr_t > & breakpoints,
                                                const std::string & sampleName )
    {
        const auto varReadBaseQualities =
            getVariantsReadBaseQualities( readPtr->getStartPos(), readPtr->getQualities(), variants, breakpoints );

        assert( variants.size() == varReadBaseQualities.size() );

        for ( std::size_t index = 0; index < variants.size(); ++index )
        {
            auto var = variants[index];
            auto qual = varReadBaseQualities[index];
            // Don't add indels at start of read. Maybe nicer code here?
            if ( var->end() > readPtr->getStartPos() )
            {
                this->addVariant( var, qual, readPtr, sampleName );

                const auto savedVar = m_variantCountsPerSample.find( var )->first;
                for ( const auto & readBreakpoint : breakpoints )
                {
                    readBreakpoint->addLocalVariant( savedVar );
                }
            }
        }

        this->addBreakpoints( breakpoints );
    }

    void VariantContainer::addVariant( varPtr_t varPtr,
                                       phred_t baseQual,
                                       io::readPtr_t readPtr,
                                       const std::string & sampleName )
    {
        auto find = m_variantCountsPerSample.find( varPtr );
        if ( find == m_variantCountsPerSample.end() )
        {
            // Add it to collection and set supporting reads as 0.
            m_variantCountsPerSample[varPtr][sampleName].m_totalVariantSupportingReads = 0;
            find = m_variantCountsPerSample.find( varPtr );
        }

        if ( baseQual >= m_minBaseQual and readPtr->getMappingQuality() >= m_minMappingQual )
        {
            ++find->second[sampleName].m_totalVariantSupportingReads;
        }

        find->first->addRead( readPtr );
    }

    //-----------------------------------------------------------------------------------------

    variantSet_t VariantContainer::getVariants() const
    {
        variantSet_t variants;
        for ( auto variantPair : m_variantCountsPerSample )
        {
            variants.insert( variantPair.first );
        }
        return variants;
    }

    void VariantContainer::addBreakpoints( const std::vector< breakpointPtr_t > & breakpoints )
    {
        for ( const auto & breakpoint : breakpoints )
        {
            if ( breakpoint->isStartBreakpoint() )
            {
                if ( m_startBreakpointLoci.find( breakpoint->pos() ) == m_startBreakpointLoci.end() )
                {
                    m_startBreakpointLoci[breakpoint->pos()] = std::make_shared< BreakpointLocus >(
                        breakpoint->contig(), breakpoint->pos(), breakpoint->isStartBreakpoint() );
                }
                m_startBreakpointLoci[breakpoint->pos()]->add( breakpoint );
            }
            else
            {
                if ( m_endBreakpointLoci.find( breakpoint->pos() ) == m_endBreakpointLoci.end() )
                {
                    m_endBreakpointLoci[breakpoint->pos()] = std::make_shared< BreakpointLocus >(
                        breakpoint->contig(), breakpoint->pos(), breakpoint->isStartBreakpoint() );
                }
                m_endBreakpointLoci[breakpoint->pos()]->add( breakpoint );
            }
        }
    }

    breakpointLocusSet_t VariantContainer::getBreakpoints() const
    {
        breakpointLocusSet_t ret;
        for ( const auto & posToBps : m_startBreakpointLoci )
        {
            if ( posToBps.second->hasMinSupport() )
            {
                ret.insert( posToBps.second );
            }
        }
        for ( const auto & posToBps : m_endBreakpointLoci )
        {
            if ( posToBps.second->hasMinSupport() )
            {
                ret.insert( posToBps.second );
            }
        }
        return ret;
    }

    //-----------------------------------------------------------------------------------------

    void VariantContainer::computeCoverage( const caller::Region & blockRegion, io::perSampleRegionsReads_t allReads )
    {
        for ( auto & variantSampleCounts : m_variantCountsPerSample )
        {
            auto varPtr = variantSampleCounts.first;
            if ( blockRegion.contains( varPtr->region() ) )
            {

                const auto overlapTester = [varPtr]( const io::Read & read )
                {
                    return varPtr->interval().overlaps( read.getMaximalReadInterval() );
                };

                const auto perSampleReadRanges = io::reduceRegionSet( allReads, varPtr->region() );

                for ( const auto & perSampleReadRange : perSampleReadRanges )
                {
                    const auto & sampleName = perSampleReadRange.first;
                    const auto & range = perSampleReadRange.second;

                    variantSampleCounts.second[sampleName].m_totalReads =
                        std::count_if( range.begin(), range.end(), overlapTester );
                }
            }
        }
    }

    //-----------------------------------------------------------------------------------------

    int64_t VariantContainer::totalReadsSupportingVariant( varPtr_t varPtr ) const
    {
        const auto variantSampleCounts = m_variantCountsPerSample.at( varPtr );
        int64_t returnValue = 0;

        for ( const auto & variantSampleCount : variantSampleCounts )
        {
            const VariantCounts & varCount = variantSampleCount.second;
            returnValue += varCount.m_totalVariantSupportingReads;
        }

        return returnValue;
    }

    //-----------------------------------------------------------------------------------------

    int64_t VariantContainer::maxReadPercentVariantCoverage( varPtr_t varPtr ) const
    {
        const auto variantSampleCounts = m_variantCountsPerSample.at( varPtr );
        int64_t maxCount = 0;

        for ( const auto & variantSampleCount : variantSampleCounts )
        {
            const VariantCounts & varCount = variantSampleCount.second;
            maxCount = std::max( maxCount, varCount.getPercentVariantCoverage() );
        }

        return maxCount;
    }
}
}
