// All content Copyright (C) 2018 Genomics plc
#include <queue>
#include "variant/breakpointVariantGenerator.hpp"
#include "variant/haplotypeGenerator.hpp"

namespace wecall
{
namespace variant
{
    BreakpointVariantGenerator::BreakpointVariantGenerator( const utils::referenceSequencePtr_t & referenceSequence,
                                                            const std::size_t kmerSize,
                                                            const std::size_t maxKmerSize,
                                                            const std::size_t minChainSupport )
        : m_defaultKmerSize( kmerSize ),
          m_maxKmerSize( maxKmerSize ),
          m_minChainSupport( minChainSupport ),
          m_refSeq( referenceSequence )
    {
    }

    assembly::SequenceGraph BreakpointVariantGenerator::processData( const caller::SetRegions & referenceRegions,
                                                                     const io::perSampleRegionsReads_t & reads ) const
    {
        for ( std::size_t kmerSize = m_defaultKmerSize; kmerSize < m_maxKmerSize; kmerSize += m_defaultKmerSize )
        {
            auto pair = this->createSequenceGraph( referenceRegions, reads, kmerSize, true );
            const auto completed = pair.second;
            if ( completed )
            {
                return pair.first;
            }
        }
        return this->createSequenceGraph( referenceRegions, reads, m_maxKmerSize, false ).first;
    }

    std::pair< assembly::SequenceGraph, bool > BreakpointVariantGenerator::createSequenceGraph(
        const caller::SetRegions & referenceRegions,
        const io::perSampleRegionsReads_t & reads,
        const std::size_t kmerSize,
        const bool disallowRepeats ) const
    {
        const std::size_t minEdgeBaseQuality = 20;
        assembly::SequenceGraph seqGraph( kmerSize, minEdgeBaseQuality );

        for ( const auto & region : referenceRegions )
        {
            const auto paddedRegion = region.getPadded( kmerSize ).getIntersect( m_refSeq->region() );
            if ( seqGraph.addReferenceSequence( m_refSeq, paddedRegion, disallowRepeats ) and disallowRepeats )
            {
                return std::make_pair( seqGraph, false );
            }
        }

        for ( const auto & readRange : reads )
        {
            for ( const auto & read : readRange.second )
            {
                const auto & readSequence = read.sequence();
                if ( readSequence.size() > seqGraph.kmerSize() )
                {
                    if ( seqGraph.addReadSequence( &readSequence, read.getQualities(), disallowRepeats ) and
                         disallowRepeats )
                    {
                        return std::make_pair( seqGraph, false );
                    }
                }
            }
        }
        return std::make_pair( seqGraph, true );
    }

    caller::Region BreakpointVariantGenerator::referenceRegionFromBreakpoint( const breakpointLocusPtr_t & bp ) const
    {
        const auto breakpointVariants = bp->getLocalVariants();

        const auto enterBpPadding = static_cast< int64_t >( m_maxKmerSize );
        const auto exitBpPadding = static_cast< int64_t >( m_defaultKmerSize );

        if ( bp->isStartLocus() )
        {
            int64_t minVarStart = bp->pos();
            for ( const auto & var : breakpointVariants )
            {
                minVarStart = std::min( minVarStart, var->start() );
            }

            int64_t refStartPos = std::max( minVarStart - enterBpPadding, m_refSeq->start() );
            int64_t redEndPos = std::min( bp->pos() + exitBpPadding, m_refSeq->end() );

            return caller::Region( bp->contig(), refStartPos, redEndPos );
        }
        else
        {
            int64_t maxVarEnd = bp->pos();
            for ( const auto & var : breakpointVariants )
            {
                maxVarEnd = std::max( maxVarEnd, var->end() );
            }

            int64_t refStartPos = std::max( bp->pos() - exitBpPadding, m_refSeq->start() );
            int64_t refEndPos = std::min( maxVarEnd + enterBpPadding, m_refSeq->end() );

            return caller::Region( bp->contig(), refStartPos, refEndPos );
        }
    }

    variantSet_t BreakpointVariantGenerator::getVariantCandidates( const breakpointLocusSet_t & breakpointLoci,
                                                                   const io::perSampleRegionsReads_t & reads ) const
    {
        caller::SetRegions referenceRegions;
        variantSet_t breakpointVariants;
        for ( const auto & bpLocus : breakpointLoci )
        {
            referenceRegions.insert( this->referenceRegionFromBreakpoint( bpLocus ) );
            const auto bpVariants = bpLocus->getLocalVariants();
            breakpointVariants.insert( bpVariants.cbegin(), bpVariants.cend() );
        }

        const auto seqGraph = this->processData( referenceRegions, reads );
        const auto chainsForVariants = seqGraph.getPathsBetweenRefNodes( m_minChainSupport );

        WECALL_LOG( DEBUG, "Assembling over " << referenceRegions );
        WECALL_LOG( DEBUG, "Created sequence graph size: " << seqGraph.kmerSize() * seqGraph.size()
                                                            << " (kmer-size * nNodes)" );
        WECALL_LOG( DEBUG, "Found nChains: " << chainsForVariants.size() );
        WECALL_LOG( SUPER_DEBUG, "Constructed seq graph:\n" << seqGraph.toString() );

        variantSet_t variants;
        for ( const auto & chain : chainsForVariants )
        {
            const auto chainVariants = chain.getVariants( m_refSeq );
            for ( auto chainVariant : chainVariants )
            {
                for ( auto rIt = breakpointVariants.rbegin(); rIt != breakpointVariants.rend(); ++rIt )
                {
                    if ( chainVariant->removable( *rIt ) )
                    {
                        chainVariant = chainVariant->remove( *rIt );
                    }
                }
                for ( auto it = breakpointVariants.begin(); it != breakpointVariants.end(); ++it )
                {
                    if ( chainVariant->removable( *it ) )
                    {
                        chainVariant = chainVariant->remove( *it );
                    }
                }

                if ( not chainVariant->empty() )
                {

                    chainVariant = chainVariant->getLeftAligned( m_refSeq->start() );
                    if ( chainVariant->isMNP() )
                    {
                        const auto split = chainVariant->split();
                        variants.insert( split.cbegin(), split.cend() );
                    }
                    else
                    {
                        variants.insert( chainVariant );
                    }
                }
            }
        }
        const auto setBreakpointVar = []( varPtr_t v )
        {
            v->setFromBreakpoint();
        };
        std::for_each( variants.cbegin(), variants.cend(), setBreakpointVar );

        return variants;
    }

    MultiLocusBreakpointVariantGenerator::MultiLocusBreakpointVariantGenerator(
        const caller::Region & region,
        const utils::referenceSequencePtr_t & referenceSequence,
        const std::size_t defaultKmerSize,
        const std::size_t maxKmerSize,
        const int largeVariantSizeDefinition )
        : m_region( region ),
          m_referenceSequence( referenceSequence ),
          m_defaultKmerSize( defaultKmerSize ),
          m_maxKmerSize( maxKmerSize ),
          m_largeVariantSizeDefinition( largeVariantSizeDefinition )
    {
    }

    std::vector< std::vector< breakpointLocusPtr_t > > BreakpointClusterer::getClusters(
        const breakpointLocusSet_t & breakpointSet ) const
    {
        breakpointLocusSet_t todo( breakpointSet.begin(), breakpointSet.end() );

        std::vector< std::vector< breakpointLocusPtr_t > > clusters;

        while ( not todo.empty() )
        {
            auto newBreakpointLocus = *todo.begin();

            todo.erase( newBreakpointLocus );

            BreakpointLocusCluster newCluster( m_minDistance );
            newCluster.push( newBreakpointLocus );

            bool notFinished = true;
            while ( notFinished )
            {
                std::vector< breakpointLocusPtr_t > related = {};
                for ( auto iterBpLocus = todo.begin(); iterBpLocus != todo.end(); ++iterBpLocus )
                {
                    if ( newCluster.isRelated( *iterBpLocus ) )
                    {
                        related.push_back( *iterBpLocus );
                    }
                }

                notFinished = ( not related.empty() );
                for ( const auto & relatedBp : related )
                {
                    todo.erase( relatedBp );
                    newCluster.push( relatedBp );
                }
            }

            clusters.push_back( newCluster.get() );
        }

        return clusters;
    }

    variantSet_t MultiLocusBreakpointVariantGenerator::getVariantCandidates(
        const breakpointLocusSet_t & breakpointSet,
        const io::perSampleRegionsReads_t & allReads ) const
    {
        const BreakpointClusterer clusterer( 30 );
        const auto clusters = clusterer.getClusters( breakpointSet );

        const std::size_t minSupport = 100;
        const std::size_t maxReadsPerSample = 5000;
        BreakpointVariantGenerator breakpointVariantGenerator( m_referenceSequence, m_defaultKmerSize, m_maxKmerSize,
                                                               minSupport );

        variantSet_t variantSet;
        // TODO(ES): More sensible method of getting kmer size.
        for ( const auto & cluster : clusters )
        {
            breakpointLocusSet_t breakpointLoci;
            caller::SetRegions breakpointReadRegions;
            for ( const auto & bpLocus : cluster )
            {
                if ( m_region.contains( bpLocus->contig(), bpLocus->pos() ) )
                {
                    breakpointLoci.insert( bpLocus );

                    const auto locusRegion = caller::Region( bpLocus->contig(), bpLocus->pos(), bpLocus->pos() )
                                                 .getPadded( m_defaultKmerSize )
                                                 .getIntersect( m_region );
                    breakpointReadRegions.insert( locusRegion );
                }
            }

            if ( not breakpointLoci.empty() )
            {
                const auto generatorReads = io::reduceRegionSet( allReads, breakpointReadRegions );

                bool tooManyReads = false;
                for ( const auto & readRange : generatorReads )
                {
                    const unsigned int nReadsThisSample =
                        std::distance( readRange.second.begin(), readRange.second.end() );
                    if ( nReadsThisSample >= maxReadsPerSample )
                    {
                        tooManyReads = true;
                    }
                }

                if ( not tooManyReads )
                {
                    const auto newVars =
                        breakpointVariantGenerator.getVariantCandidates( breakpointLoci, generatorReads );
                    for ( const auto & newVar : newVars )
                    {
                        if ( newVar->sequenceLength() >= m_largeVariantSizeDefinition or
                             newVar->sequenceLengthInRef() >= m_largeVariantSizeDefinition )
                        {
                            variantSet.insert( newVar );
                        }
                    }
                }
                else
                {
                    WECALL_LOG( WARNING, "Skipping assembly on " << breakpointReadRegions
                                                                  << " due to too many reads" );
                }
            }
        }

        return variantSet;
    }
}
}
