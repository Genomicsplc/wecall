// All content Copyright (C) 2018 Genomics plc
#include "caller/job.hpp"
#include "caller/params.hpp"
#include "caller/annotation.hpp"
#include "caller/diploid/diploid.hpp"
#include "caller/regionUtils.hpp"
#include "caller/diploid/referenceCalling.hpp"
#include "caller/alignPhasing.hpp"

#include "variant/variantGenerator.hpp"
#include "variant/variantFilter.hpp"
#include "variant/variantContainer.hpp"
#include "variant/type/variant.hpp"
#include "variant/clustering.hpp"
#include "variant/haplotypeGenerator.hpp"
#include "variant/type/variant.hpp"
#include "variant/breakpointVariantGenerator.hpp"

#include "io/bamFile.hpp"
#include "io/readRange.hpp"
#include "vcf/reader.hpp"

#include "readrecalibration/readRecalibration.hpp"

#include "utils/exceptions.hpp"
#include "utils/logging.hpp"
#include "utils/flatten.hpp"

#include <algorithm>
#include <limits>
#include <boost/filesystem/operations.hpp>

using namespace echidna::io;
using echidna::variant::varPtr_t;

namespace echidna
{
namespace caller
{
    //-----------------------------------------------------------------------------------------

    Job::Job( caller::params::Application applicationParams,
              caller::params::Data dataParams,
              caller::params::System systemParams,
              caller::params::PrivateSystem privateSystemParams,
              caller::params::Filters filterParams,
              caller::params::PrivateCalling privateCallingParams,
              caller::params::Calling callingParams,
              caller::params::PrivateData privateDataParams )
        : m_applicationParams( applicationParams ),
          m_dataParams( dataParams ),
          m_systemParams( systemParams ),
          m_privateSystemParams( privateSystemParams ),
          m_filterParams( filterParams ),
          m_privateCallingParams( privateCallingParams ),
          m_callingParams( callingParams ),
          m_privateDataParams( privateDataParams ),
          m_candidateVCFFile(
              privateDataParams.m_candidateVariantsFile.empty()
                  ? nullptr
                  : std::make_shared< io::TabixVCFFile >( privateDataParams.m_candidateVariantsFile,
                                                          privateDataParams.m_candidateVariantsFile + ".tbi" ) ),
          m_genotypeAllelesVCFFile(
              privateDataParams.genotypingMode()
                  ? std::make_shared< io::TabixVCFFile >( privateDataParams.genotypeAllelesFile(),
                                                          privateDataParams.genotypeAllelesFile() + ".tbi" )
                  : nullptr ),
          m_intermediateOutputWriter( dataParams.inputDataSources(), privateDataParams.m_intermediateRecalibFileStem ),
          m_variantSoftFilterBank( privateCallingParams.m_varFilterIDs,
                                   callingParams.m_minAlleleBiasP,
                                   callingParams.m_minStrandBiasP,
                                   callingParams.m_minAllelePlusStrandBiasP,
                                   callingParams.m_minRMSMappingQ,
                                   callingParams.m_minSNPQOverDepth,
                                   callingParams.m_minIndelQOverDepth,
                                   callingParams.m_minBadReadsScore,
                                   callingParams.m_minCallQual ),
          m_readDataReader( systemParams,
                            privateSystemParams,
                            filterParams,
                            privateSystemParams.m_biteSize,
                            dataParams.inputDataSources() ),
          m_vcOut( dataParams.outputDataSink(), dataParams.outputRefCalls(), callingParams.m_outputPhasedGenotypes ),
          m_ref( dataParams.refFile() ),
          // TODO(ES): Tie together contig, calling and output regions together into nice container.
          m_outputRegions( utils::functional::flatten( dataParams.dataRegions() ) ),
          m_callingRegions( utils::functional::flatten(
              caller::padPartitionedRegions( dataParams.dataRegions(), m_privateCallingParams.m_regionPadding ) ) ),
          m_model( callingParams.m_badReadsWindowSize,
                   privateCallingParams.m_maxHaplotypesPerCluster,
                   m_readDataReader.getSampleNames() ),
          m_variantCallBuilder( dataParams.outputRefCalls(),
                                callingParams.m_outputPhasedGenotypes,
                                privateCallingParams.m_allVariants,
                                privateDataParams.genotypingMode(),
                                m_readDataReader.getSampleNames().size() )
    {
    }

    //-----------------------------------------------------------------------------------------

    void Job::process()
    {
        ECHIDNA_LOG( DEBUG, "Writing VCF header" );

        std::vector< vcf::FilterDesc > filterDescs = m_variantSoftFilterBank.getFilterDescs();
        filterDescs.emplace_back( vcf::filter::NC_key,
                                  "Not called: Indicates a variant that was not positively genotyped in any sample." );

        m_vcOut.writeHeader( m_dataParams.outputFormat(), m_applicationParams, m_dataParams.refFile(),
                             m_readDataReader.getSampleNames(), filterDescs,
                             caller::getContigsFromRegions( m_outputRegions, m_ref.indexFile().contigs() ) );

        // Currently defaulting to same ploidy all samples.
        const std::vector< std::size_t > ploidyPerSample( m_readDataReader.getSampleNames().size(),
                                                          m_privateCallingParams.m_ploidy );

        for ( const Region & region : m_callingRegions )
        {
            const auto contig = region.contig();
            const auto refContigs = m_ref.indexFile().contigs();
            const auto contigLookup = refContigs.find( contig );

            if ( contigLookup == refContigs.end() )
            {
                continue;
            }

            const auto contigRegion = Region( contig, contigLookup->second );

            if ( not contigRegion.overlaps( region ) )
            {
                continue;
            }

            const auto callingRegion = contigRegion.getIntersect( region );

            if ( callingRegion.size() == 0 )
            {
                continue;
            }

            m_vcOut.contig( contig );

            // Get required reference.
            const auto maxReadLength =
                1000;  // TODO made sure Read construction doesn't fail it it goes byond this limit
            const auto assemblePadding = 3 * params::defaults::maxBreakpointKmerSize;
            const caller::Region paddedRefRegion = callingRegion.getPadded( maxReadLength + assemblePadding );
            m_ref.cacheSequence( paddedRefRegion );
            const auto referenceSequence =
                std::make_shared< utils::ReferenceSequence >( m_ref.getSequence( paddedRefRegion ) );

            auto blockIterator = m_readDataReader.readRegion( callingRegion, referenceSequence );

            while ( auto readDataset = blockIterator.getReadDatasetForNextBlock() )
            {
                const auto block = readDataset->region();

                if ( not readDataset->isEmpty() or m_dataParams.outputRefCalls() or
                     m_privateDataParams.genotypingMode() )
                {
                    const auto actualEnd = processBlock( readDataset, blockIterator.isLastBlock(), ploidyPerSample );

                    if ( actualEnd < block.end() )
                    {
                        blockIterator.chopCurrentBlock( actualEnd );
                    }
                }
                else
                {
                    ECHIDNA_LOG( INFO, "Skipping " << block << " due to no read data." );
                }
            }
        }

        ECHIDNA_LOG( INFO, "Job completed successfully." );
    }

    //-----------------------------------------------------------------------------------------

    std::vector< variant::VariantCluster > Job::generateVariantClustersInBlock(
        const caller::Region & blockRegion,
        io::perSampleRegionsReads_t allReads,
        utils::referenceSequencePtr_t referenceSequence )
    {
        const variant::VariantGenerator varGen( referenceSequence, m_privateCallingParams.m_minBaseQual,
                                                m_filterParams.m_readMappingFilterQ );

        auto varContainer = varGen.generateVariantsFromReads( allReads );

        if ( m_candidateVCFFile != nullptr )
        {
            for ( const auto & record : m_candidateVCFFile->fetch( blockRegion ) )
            {
                const auto variants = record.getVariants( referenceSequence );
                const auto priors = caller::getPriorsFromInfo( record.m_info );
                for ( std::size_t index = 0; index < variants.size(); ++index )
                {
                    varContainer.addCandidateVariant( variants[index], priors[index] );
                }
            }
        }

        if ( m_genotypeAllelesVCFFile != nullptr )
        {
            for ( const auto & record : m_genotypeAllelesVCFFile->fetch( blockRegion ) )
            {
                const auto variants = record.getVariants( referenceSequence );
                for ( std::size_t index = 0; index < variants.size(); ++index )
                {
                    varContainer.addGenotypingVariant( variants[index] );
                }
            }
        }

        varContainer.computeCoverage( blockRegion, allReads );

        const auto varFilter = variant::VariantFilter( m_privateCallingParams.m_minReadsPerVar,
                                                       m_privateCallingParams.m_perSamPercentReadsPerVar );

        auto variants = varFilter.getSortedFilteredVariants( blockRegion, varContainer );

        if ( m_callingParams.m_turnOnLargeVariantCalls )
        {
            const auto breakpoints = varContainer.getBreakpoints();

            const variant::MultiLocusBreakpointVariantGenerator multiLocusBreakpointVariantGenerator(
                blockRegion, referenceSequence, params::defaults::defaultBreakpointKmerSize,
                params::defaults::maxBreakpointKmerSize, m_privateCallingParams.m_largeVariantSizeDefinition );
            const auto bpVariants = multiLocusBreakpointVariantGenerator.getVariantCandidates( breakpoints, allReads );

            for ( const auto & bpVar : bpVariants )
            {
                // This is a problem that small CNVs know no boundaries...
                if ( blockRegion.contains( bpVar->region() ) )
                {
                    variants.insert( bpVar );
                }
            }
        }

        return variant::generateMergedClusters( variants, m_privateCallingParams.m_minClusterDist,
                                                m_privateCallingParams.m_maxClusterDist,
                                                m_privateCallingParams.m_maxClusterSize, blockRegion, referenceSequence,
                                                m_privateCallingParams.m_maxClusterVariantCombinations,
                                                m_privateCallingParams.m_minReadsToMakeCombinationClaim );
    }

    //-----------------------------------------------------------------------------------------

    int64_t Job::processBlock( io::readDataset_t readDataset,
                               bool lastBlock,
                               const std::vector< std::size_t > & ploidyPerSample )
    {
        auto blockRegion = readDataset->region();
        ECHIDNA_LOG( INFO, "Processing:\t" << blockRegion );
        // TODO: Calibrate min value below.
        const auto allReads = readDataset->getAllReads( m_filterParams.m_readMappingFilterQ );

        // Get required reference.
        const auto maxReadLength = io::perSampleMaxAlignedReadLength( allReads );
        const auto assemblePadding = 3 * params::defaults::maxBreakpointKmerSize;
        const caller::Region paddedRefRegion = blockRegion.getPadded( maxReadLength + assemblePadding );
        m_ref.cacheSequence( paddedRefRegion );

        const auto referenceSequence =
            std::make_shared< utils::ReferenceSequence >( m_ref.getSequence( paddedRefRegion ) );

        auto clusters = this->generateVariantClustersInBlock( blockRegion, allReads, referenceSequence );
        if ( not lastBlock )
        {
            variant::removeClustersNearBlockEnd( blockRegion, clusters, m_privateCallingParams.m_maxClusterDist );
        }

        ECHIDNA_LOG( INFO, "Calling variants:\t" << blockRegion );

        if ( clusters.empty() )
        {
            this->callReference( readDataset->region().contig(), blockRegion.start(), blockRegion.end(), readDataset,
                                 ploidyPerSample );
        }
        else
        {
            // init with first cluster
            std::vector< variant::VariantCluster >::iterator clusterIterator = clusters.begin();
            auto firstCluster = *clusterIterator;
            auto callsPrevCluster =
                this->processBigCluster( firstCluster, readDataset, referenceSequence, ploidyPerSample );
            this->writeCallsForCluster( callsPrevCluster, firstCluster, blockRegion.start(), readDataset,
                                        ploidyPerSample );
            Region regionPrevCluster = firstCluster.region();

            // iterate over remaining clusters
            clusterIterator++;
            for ( ; clusterIterator != clusters.cend(); ++clusterIterator )
            {
                // call variants in cluster
                const auto cluster = *clusterIterator;
                auto calls = this->processBigCluster( cluster, readDataset, referenceSequence, ploidyPerSample );

                callVector_t aligned_calls;

                // align phasing
                if ( callsPrevCluster.empty() )
                {
                    aligned_calls = calls;
                }
                else
                {
                    // get padded reference sequence for phasing
                    const Region combinedRegion =
                        Region( cluster.region().contig(), regionPrevCluster.start(), cluster.region().end() );
                    const io::perSampleRegionsReads_t combinedRegionReads =
                        readDataset->getRegionsReads( combinedRegion, m_filterParams.m_readMappingFilterQ );
                    const variant::VariantCluster combinedCluster = variant::VariantCluster( {}, combinedRegion );
                    const auto paddedReferenceSequence =
                        this->getReferenceForCluster( combinedCluster, combinedRegionReads, referenceSequence );

                    // attempt to phase clusters
                    alignPhasingBetweenClusters( calls, callsPrevCluster, cluster.region(), regionPrevCluster,
                                                 combinedRegion, paddedReferenceSequence, ploidyPerSample,
                                                 combinedRegionReads, readDataset->getSampleNames() );
                    aligned_calls = calls;
                }

                this->writeCallsForCluster( aligned_calls, cluster, regionPrevCluster.end(), readDataset,
                                            ploidyPerSample );

                // move one cluster forward
                callsPrevCluster = calls;
                regionPrevCluster = cluster.region();
            }

            // call block at the end
            this->callReference( readDataset->region().contig(), regionPrevCluster.end(), blockRegion.end(),
                                 readDataset, ploidyPerSample );
        }
        return blockRegion.end();
    }

    //-----------------------------------------------------------------------------------------

    void Job::writeCallsForCluster( const callVector_t calls,
                                    const variant::VariantCluster cluster,
                                    const int64_t refStart,
                                    io::readDataset_t readDataset,
                                    const std::vector< std::size_t > & ploidyPerSample )
    {
        callVector_t outputCalls = this->filterOutputCalls( readDataset->region().contig(), calls );
        // output reference calls between previous and current cluster
        this->callReference( readDataset->region().contig(), refStart, cluster.zeroIndexedVcfStart( outputCalls ),
                             readDataset, ploidyPerSample );

        // write calls (variant + ref calls) for current cluster
        m_variantSoftFilterBank.applyFilterAnnotation( outputCalls );
        m_vcOut.writeCallSet( m_ref, outputCalls );
    }

    //-----------------------------------------------------------------------------------------

    void Job::callReference( const std::string & contig,
                             int64_t start,
                             int64_t end,
                             io::readDataset_t readDataset,
                             const std::vector< std::size_t > & ploidyPerSample )
    {
        if ( m_dataParams.outputRefCalls() and start < end )
        {
            const caller::Region refInterval( contig, start, end );
            const auto reads = readDataset->getRegionsReads( refInterval, m_filterParams.m_readMappingFilterQ );

            callVector_t calls;
            const auto maxUncalledVaraintQuality = 0.0;
            const int64_t maxRefCallSize = m_dataParams.maxRefCallSize();

            for ( auto begBlock = refInterval.start(); begBlock < refInterval.end(); begBlock += maxRefCallSize )
            {
                const int64_t endBlock = std::min( begBlock + maxRefCallSize, refInterval.end() );
                const caller::Region region( contig, utils::Interval( begBlock, endBlock ) );
                const auto chunkedRefCalls =
                    caller::model::buildRefCall( region, reads, maxUncalledVaraintQuality, ploidyPerSample,
                                                 m_privateCallingParams.m_referenceCallQualityDeltaThreshold );
                calls.insert( calls.end(), chunkedRefCalls.begin(), chunkedRefCalls.end() );
            }

            const auto outputCalls = this->filterOutputCalls( contig, calls );

            m_vcOut.writeCallSet( m_ref, outputCalls );
        }
    }

    //-----------------------------------------------------------------------------------------

    void Job::recalibrateBaseQualities( const io::perSampleRegionsReads_t & readRangesPerSample,
                                        const caller::Region & region ) const
    {
        const auto errorCorrectionParameters = corrector::ErrorCorrectionParameters();
        corrector::recalibrateReads( readRangesPerSample, m_ref, region, errorCorrectionParameters );
        m_intermediateOutputWriter.writeReads( readRangesPerSample, region.contig() );
    }

    std::vector< varPtr_t > getCandidateVariants( const std::vector< varPtr_t > & variants,
                                                  const variant::variantSet_t & newVariants )
    {
        if ( newVariants.size() == 0 )
        {
            return variants;
        }
        else
        {
            variant::variantSet_t variantsWithMNPs( variants.cbegin(), variants.cend() );
            variantsWithMNPs.insert( newVariants.cbegin(), newVariants.cend() );
            return std::vector< varPtr_t >( variantsWithMNPs.cbegin(), variantsWithMNPs.cend() );
        }
    }

    std::vector< std::size_t > getReducedPloidies( const callVector_t & calls,
                                                   const std::vector< std::size_t > & defaultPloidies,
                                                   const caller::Region & thisRegion )
    {
        auto reducedPloidies = defaultPloidies;
        for ( std::size_t sampleIndex = 0; sampleIndex < defaultPloidies.size(); ++sampleIndex )
        {
            auto & currentPloidy = reducedPloidies[sampleIndex];
            for ( const auto & call : calls )
            {
                if ( not call.isRefCall() and call.var->region().contains( thisRegion ) )
                {
                    const auto & genotypes = call.samples[sampleIndex].genotypeCalls;
                    for ( const auto & genotype : genotypes )
                    {
                        if ( genotype == caller::Call::VAR )
                        {
                            --currentPloidy;
                        }
                    }
                }
            }
        }

        return reducedPloidies;
    }

    callVector_t Job::processBigCluster( const variant::VariantCluster & cluster,
                                         io::readDataset_t readDataset,
                                         const utils::referenceSequencePtr_t & blockReferenceSequence,
                                         const std::vector< std::size_t > & ploidyPerSample )
    {
        variant::setDefaultPriors( cluster.variants() );

        const auto bigClusterReads =
            readDataset->getRegionsReads( cluster.region(), m_filterParams.m_readMappingFilterQ );
        const auto allReads = readDataset->getAllReads( m_filterParams.m_readMappingFilterQ );

        bool hasLargeVariant = false;
        const auto largeVariantThreshold = m_privateCallingParams.m_largeVariantClusterThreshold;
        if ( m_callingParams.m_turnOnLargeVariantCalls )
        {
            hasLargeVariant = std::any_of( cluster.variants().cbegin(), cluster.variants().cend(),
                                           [largeVariantThreshold]( const varPtr_t & var )
                                           {
                                               return var->isLargeVariant( largeVariantThreshold );
                                           } );
        }

        if ( not hasLargeVariant )
        {
            return processCluster( cluster, bigClusterReads, allReads, blockReferenceSequence, ploidyPerSample );
        }
        else
        {
            const auto data = cluster.buildSubClusters( m_privateCallingParams.m_maxClusterDist );
            auto lvcCluster = std::get< 0 >( data );
            auto smallVariantClustersNotTouchingLargeVariants = std::get< 1 >( data );

            lvcCluster.computeVariantCombinations(
                m_privateCallingParams.m_minReadsToMakeCombinationClaim, m_privateCallingParams.m_maxClusterDist,
                m_privateCallingParams.m_maxClusterVariantCombinations, blockReferenceSequence );
            for ( auto & miniCluster : smallVariantClustersNotTouchingLargeVariants )
            {
                miniCluster.computeVariantCombinations(
                    m_privateCallingParams.m_minReadsToMakeCombinationClaim, m_privateCallingParams.m_maxClusterDist,
                    m_privateCallingParams.m_maxClusterVariantCombinations, blockReferenceSequence );
            }

            const auto lvcReads = io::reduceRegionSet( bigClusterReads, lvcCluster.readRegions() );
            auto lvcCalls = processCluster( lvcCluster, lvcReads, allReads, blockReferenceSequence, ploidyPerSample );

            callVector_t allCalls = lvcCalls;
            for ( const auto & leftOverCluster : smallVariantClustersNotTouchingLargeVariants )
            {
                const auto clusterReads = io::reduceRegionSet( bigClusterReads, leftOverCluster.readRegions() );
                const auto thisAreasPloidy = getReducedPloidies( lvcCalls, ploidyPerSample, leftOverCluster.region() );
                auto leftOverCalls =
                    processCluster( leftOverCluster, clusterReads, allReads, blockReferenceSequence, thisAreasPloidy );
                const auto callMerger = lvcMerger();
                auto leftOverCallsWithoutReference =
                    callMerger.removeReferenceCalls( leftOverCalls, thisAreasPloidy, ploidyPerSample );
                callMerger.mergeAndCorrectGenotypes( leftOverCallsWithoutReference, lvcCalls, thisAreasPloidy,
                                                     ploidyPerSample );
                allCalls.insert( allCalls.end(), leftOverCallsWithoutReference.begin(),
                                 leftOverCallsWithoutReference.end() );
            }

            std::sort( allCalls.begin(), allCalls.end(), caller::CallComp() );

            return allCalls;
        }
    }

    utils::referenceSequencePtr_t Job::getReferenceForCluster(
        const variant::VariantCluster & cluster,
        const io::perSampleRegionsReads_t & regionReads,
        const utils::referenceSequencePtr_t & blockReferenceSequence ) const
    {
        const int64_t minRefPadding = 40L;  // TODO: Handle this more safely. Needed for kmer building later in code.
        const auto maxReadLength = io::perSampleMaxReadCigarLength( regionReads );
        const auto refPaddingSize = std::max( minRefPadding, maxReadLength + constants::needlemanWunschPadding );

        return std::make_shared< utils::ReferenceSequence >(
            blockReferenceSequence->getPaddedSequence( cluster.region(), cluster.paddedRegion(), refPaddingSize ) );
    }

    //-----------------------------------------------------------------------------------------

    callVector_t Job::processCluster( const variant::VariantCluster & cluster,
                                      const io::perSampleRegionsReads_t & regionReads,
                                      const io::perSampleRegionsReads_t & allReads,
                                      const utils::referenceSequencePtr_t & blockReferenceSequence,
                                      const std::vector< std::size_t > & ploidyPerSample )
    {
        ECHIDNA_LOG( DEBUG, "Processing: " << cluster.toString() );

        if ( m_callingParams.m_recalibrateBaseQs )
        {
            this->recalibrateBaseQualities( regionReads, cluster.region() );
        }

        const auto paddedRefSequence = this->getReferenceForCluster( cluster, regionReads, blockReferenceSequence );

        variant::HaplotypeVector haplotypes( cluster.readRegions(), paddedRefSequence );

        if ( cluster.allCombinationsComputed() )
        {
            for ( const auto & varCombo : cluster.variantCombinations() )
            {
                haplotypes.push_back( variant::variantSet_t( varCombo.cbegin(), varCombo.cend() ) );
            }
            haplotypes.sort();
            haplotypes.merge();
        }
        else
        {
            ECHIDNA_LOG( DEBUG, "Re-clustering " << cluster.region()
                                                 << " with nVariants: " << cluster.variants().size() );
            const variant::AlignmentHaplotypeGenerator hapGen(
                cluster.variants(), cluster.readRegions(), regionReads, paddedRefSequence,
                m_privateCallingParams.m_maxHaplotypesPerCluster,
                m_privateCallingParams.m_minReadsToMakeCombinationClaim );

            haplotypes = hapGen.generateHaplotypes();
            if ( haplotypes.size() <= 1 )
            {
                ECHIDNA_LOG( WARNING, "Skipping region: " << cluster.region() << " with " << cluster.variants().size()
                                                          << " variants." );
                return {};
            }
        }
        ECHIDNA_LOG( DEBUG, "Processing " << haplotypes.size() << " haplotypes" );

        std::vector< varPtr_t > candidateVariants;
        if ( m_privateCallingParams.m_normalizeVariantCalls )
        {
            const auto variants = haplotypes.withNormalizedVariants();
            candidateVariants = std::vector< varPtr_t >( variants.cbegin(), variants.cend() );
            variant::setDefaultPriors( candidateVariants );
        }
        else if ( m_callingParams.m_allowMNPCalls )
        {
            const auto newVariants = haplotypes.withMNPs();
            candidateVariants = getCandidateVariants( cluster.variants(), newVariants );

            variant::setDefaultPriors( std::vector< varPtr_t >( newVariants.cbegin(), newVariants.cend() ) );
        }
        else
        {
            candidateVariants = cluster.variants();
        }

        const auto results = m_model.getResults( regionReads, haplotypes, candidateVariants, ploidyPerSample );

        const auto calls = m_variantCallBuilder.getAnnotatedVariantCalls(
            cluster.region().start(), haplotypes.region(), regionReads, allReads, candidateVariants,
            results.getVariantQualities(), results.getVariantMetadata(), results.getCalledGenotypes(),
            results.getGenotypeMetadata(), ploidyPerSample, haplotypes,
            m_privateCallingParams.m_referenceCallQualityDeltaThreshold );

        return calls;
    }

    callVector_t Job::filterOutputCalls( const std::string & contig, const callVector_t & calls ) const
    {
        callVector_t outputCalls;
        for ( const auto & call : calls )
        {
            const Region callRegion( contig, call.interval );  // One bp before variant region for indels.
            const auto overlapCall = [callRegion]( const Region & outputRegion )
            {
                return outputRegion.overlaps( callRegion );
            };

            // TODO: Replace with better handling than linear scan.
            if ( std::any_of( m_outputRegions.begin(), m_outputRegions.end(), overlapCall ) )
            {
                outputCalls.push_back( call );
            }
        }

        return outputCalls;
    }
}
}
