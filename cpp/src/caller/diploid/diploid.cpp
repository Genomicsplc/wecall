// All content Copyright (C) 2018 Genomics plc
#include <algorithm>

#include "caller/diploid/diploid.hpp"
#include "caller/diploid/diploidAnnotate.hpp"
#include "caller/diploid/variantQualityCalculator.hpp"
#include "caller/diploid/genotypeUtils.hpp"
#include "caller/diploid/referenceCalling.hpp"
#include "variant/genotype.hpp"
#include "mapping/hashMapper.hpp"
#include "alignment/aligner.hpp"
#include "variant/genotype.hpp"
#include "alignment/galign.hpp"
#include "stats/functions.hpp"
#include "stats/models.hpp"
#include "io/readDataSet.hpp"
#include "io/readSummaries.hpp"
#include "utils/logging.hpp"
#include "utils/exceptions.hpp"
#include "utils/bestScoreSelector.hpp"
#include "utils/median.hpp"
#include "io/readUtils.hpp"
#include "caller/haplotypeLikelihoods.hpp"
#include "caller/metadata.hpp"

namespace echidna
{
namespace caller
{
    namespace model
    {

        genoCalls_t getGenoCallsForVar( const variant::HaplotypeVector & haplotypes,
                                        const std::vector< variant::genotypePtr_t > & calledGenotypes,
                                        const variant::varPtr_t & var,
                                        const std::size_t nSamples,
                                        const std::vector< std::size_t > & ploidy )
        {
            genoCalls_t genoCalls( nSamples );

            for ( std::size_t sampleIndex = 0; sampleIndex < nSamples; ++sampleIndex )
            {
                const auto genotype = calledGenotypes[sampleIndex];

                if ( genotype == nullptr )
                {
                    genoCalls[sampleIndex] = genoCall_t( ploidy[sampleIndex], caller::Call::UNKNOWN );
                }
                else
                {
                    genoCalls[sampleIndex] = genotype->getGenotypeCall( haplotypes, var );
                }
            }

            return genoCalls;
        }

        VCFCallVectorBuilder::VCFCallVectorBuilder( const bool makeRefCalls,
                                                    const bool outputPhasedCalls,
                                                    const bool outputAllVariants,
                                                    const bool genotypingMode,
                                                    const std::size_t nSamples )
            : m_makeRefCalls( makeRefCalls ),
              m_outputPhasedCalls( outputPhasedCalls ),
              m_outputAllVariants( outputAllVariants ),
              m_genotypingMode( genotypingMode ),
              m_nSamples( nSamples )
        {
        }

        callVector_t VCFCallVectorBuilder::getAnnotatedVariantCalls(
            const int64_t refStartPos,
            const caller::Region & region,
            const io::perSampleRegionsReads_t & readRangesPerSample,
            const io::perSampleRegionsReads_t & allReads,
            const std::vector< variant::varPtr_t > & variants,
            const std::vector< double > & variantQualities,
            const std::vector< std::vector< VariantMetadata > > & variantAnnotationPerSample,
            const std::vector< variant::genotypePtr_t > & calledGenotypes,
            const std::vector< GenotypeMetadata > & genotypeMetadataPerSample,
            const std::vector< std::size_t > & ploidyPerSample,
            const variant::HaplotypeVector & haplotypes,
            const double referenceCallQualityDeltaThreshold ) const
        {
            callVector_t calls;
            double maxUncalledVarQ = 0.0;
            const auto phaseSetID = region.start();
            int64_t currentRefPos = refStartPos;
            const int64_t refEndPosition = region.end();

            for ( std::size_t varIndex = 0; varIndex < variants.size(); ++varIndex )
            {
                const auto & var = variants[varIndex];
                const auto variantQuality = variantQualities[varIndex];
                const auto genoCalls =
                    getGenoCallsForVar( haplotypes, calledGenotypes, var, m_nSamples, ploidyPerSample );
                const auto variantCalled =
                    caller::variantCalled( genoCalls ) and variantQuality >= constants::minAllowedQualityScore;

                bool outputVariant = m_genotypingMode ? var->isGenotypingVariant() : variantCalled;

                if ( m_outputAllVariants or outputVariant )
                {
                    if ( m_makeRefCalls and var->zeroIndexedVcfPosition() > currentRefPos )
                    {
                        const caller::Region refRegion(
                            region.contig(), utils::Interval( currentRefPos, var->zeroIndexedVcfPosition() ) );
                        const auto chunkedRefCalls = caller::model::buildRefCall(
                            refRegion, allReads, maxUncalledVarQ, ploidyPerSample, referenceCallQualityDeltaThreshold );
                        calls.insert( calls.end(), chunkedRefCalls.begin(), chunkedRefCalls.end() );
                        maxUncalledVarQ = 0.0;  // Reset
                    }

                    // Create, annotate and add a variant call to the set
                    const utils::Interval vcfInterval( var->zeroIndexedVcfPosition(), var->end() );
                    Call call( var, vcfInterval, variantQuality, m_nSamples, genoCalls );

                    caller::annotate::annotate( call, varIndex, variantAnnotationPerSample, genotypeMetadataPerSample,
                                                m_outputPhasedCalls, phaseSetID );

                    if ( not variantCalled )
                    {
                        call.filters.insert( "NC" );
                    }

                    calls.push_back( call );

                    // Move the reference position on past this variant
                    currentRefPos = std::max( var->end(), currentRefPos );
                }
                else
                {
                    if ( m_makeRefCalls )
                    {
                        maxUncalledVarQ = std::max( maxUncalledVarQ, variantQuality );
                    }
                }
            }

            // May need a final REFCALL at the end of the set

            if ( m_makeRefCalls and refEndPosition > currentRefPos )
            {
                const caller::Region refRegion( region.contig(), utils::Interval( currentRefPos, refEndPosition ) );
                const auto chunkedRefCalls = caller::model::buildRefCall(
                    refRegion, allReads, maxUncalledVarQ, ploidyPerSample, referenceCallQualityDeltaThreshold );
                calls.insert( calls.end(), chunkedRefCalls.begin(), chunkedRefCalls.end() );
            }

            ECHIDNA_LOG( DEBUG, "Called " << calls.size() << " variants" );
            return calls;
        }

        Model::Model( const int badReadsWindowSize,
                      const std::size_t maxHaplotypesPerCluster,
                      const std::vector< std::string > & samples )
            : m_badReadsWindowSize( badReadsWindowSize ),
              m_maxHaplotypesPerCluster( maxHaplotypesPerCluster ),
              m_samples( samples )
        {
        }

        void Model::computeResultsPerSample( const std::string & sample,
                                             const std::size_t ploidy,
                                             const io::RegionsReads & readRange,
                                             const std::vector< variant::varPtr_t > & variants,
                                             const variant::HaplotypeVector & mergedHaplotypes,
                                             std::map< std::string, variant::GenotypeVector > & perSampleGenotypes,
                                             variant::genotypePtr_t & calledGenotype,
                                             GenotypeMetadata & genotypeMetadata,
                                             std::vector< double > & genotypeLikelihoods,
                                             std::vector< double > & totalHaplotypeFrequencies,
                                             std::vector< VariantMetadata > & variantAnnotation ) const
        {
            const auto haplotypeLikelihoods = computeHaplotypeLikelihoods( mergedHaplotypes, readRange );
            const bool hasReadData = haplotypeLikelihoods.size1() > 0;
            const auto haplotypeFrequencies = caller::computeHaplotypeFrequencies( haplotypeLikelihoods );

            if ( false )
            {
                // print haplotype frequencies (frequency of haplotype in population estimated from read data)
                ECHIDNA_LOG( SUPER_DEBUG, "Haplotype frequencies for " + sample + ":-" );
                for ( std::size_t haplotypeIndex = 0; haplotypeIndex < haplotypeFrequencies.size(); ++haplotypeIndex )
                {
                    ECHIDNA_LOG( SUPER_DEBUG, std::to_string( haplotypeFrequencies[haplotypeIndex] ) + " for " +
                                                  mergedHaplotypes[haplotypeIndex].toString() );
                }
            }

            for ( std::size_t haplotypeIndex = 0; haplotypeIndex != mergedHaplotypes.size(); ++haplotypeIndex )
            {
                totalHaplotypeFrequencies[haplotypeIndex] += haplotypeFrequencies[haplotypeIndex];
            }

            variantAnnotation = annotate::computeVariantAnnotation( haplotypeLikelihoods, readRange,
                                                                    m_badReadsWindowSize, variants, mergedHaplotypes );

            const auto bestIndicies = utils::indiciesWithHighestValues< double >(
                haplotypeFrequencies, this->m_maxHaplotypesPerCluster, 1.0, 1.0e5 );

            perSampleGenotypes.emplace( std::piecewise_construct, std::forward_as_tuple( sample ),
                                        std::forward_as_tuple( ploidy, mergedHaplotypes, bestIndicies, variants ) );

            const auto & genotypes = perSampleGenotypes.at( sample );

            if ( hasReadData and genotypes.size() > 0 )
            {
                genotypeLikelihoods = computeGenotypeLikelihoods( genotypes, haplotypeLikelihoods, mergedHaplotypes );

                const auto indexOfMax = utils::indexOfHighestValue( genotypeLikelihoods );

                calledGenotype = genotypes[indexOfMax];

                genotypeMetadata = annotate::computeGenotypeMetaData( indexOfMax, genotypes, genotypeLikelihoods );

                ECHIDNA_LOG( SUPER_DEBUG, "Called genotype for " << sample << ":- "
                                                                 << calledGenotype->toString( mergedHaplotypes ) );
            }
            else
            {
                genotypeLikelihoods = std::vector< double >( genotypes.size(), constants::unknownValue );
            }

            for ( std::size_t variantIndex = 0; variantIndex < variants.size(); ++variantIndex )
            {
                if ( hasReadData and genotypes.size() > 0 )
                {
                    variantAnnotation[variantIndex].genotypeLikelihoods =
                        annotate::get_RR_RA_AA_Likelihoods_as_phred_scores( variants[variantIndex], mergedHaplotypes,
                                                                            genotypes, genotypeLikelihoods );
                }
                else
                {
                    variantAnnotation[variantIndex].genotypeLikelihoods = {constants::unknownValue};
                }
            }
        }

        //-----------------------------------------------------------------------------------------

        ModelResults Model::getResults( const io::perSampleRegionsReads_t & readRangesPerSample,
                                        const variant::HaplotypeVector & mergedHaplotypes,
                                        const std::vector< variant::varPtr_t > & variants,
                                        const std::vector< std::size_t > & perSamplePloidy ) const
        {
            std::vector< std::vector< double > > genotypeLikelihoodsAllSamples( m_samples.size() );
            std::vector< variant::genotypePtr_t > calledGenotypes( m_samples.size(), nullptr );
            std::vector< std::vector< VariantMetadata > > variantAnnotationPerSample( m_samples.size() );
            std::vector< GenotypeMetadata > genotypeMetadataPerSample( m_samples.size(), GenotypeMetadata() );
            std::vector< double > totalHaplotypeFrequencies( mergedHaplotypes.size(), 0 );

            std::map< std::string, variant::GenotypeVector > perSampleGenotypes;

            for ( std::size_t sampleIndex = 0; sampleIndex < m_samples.size(); ++sampleIndex )
            {
                const auto & sample = m_samples[sampleIndex];
                const auto & readRange = readRangesPerSample.at( sample );
                const auto ploidy = perSamplePloidy[sampleIndex];
                this->computeResultsPerSample( sample, ploidy, readRange, variants, mergedHaplotypes,
                                               perSampleGenotypes, calledGenotypes[sampleIndex],
                                               genotypeMetadataPerSample[sampleIndex],
                                               genotypeLikelihoodsAllSamples[sampleIndex], totalHaplotypeFrequencies,
                                               variantAnnotationPerSample[sampleIndex] );
            }

            const double haplotypeFrequencySum =
                std::accumulate( totalHaplotypeFrequencies.cbegin(), totalHaplotypeFrequencies.cend(), 0.0 );

            for ( auto && frequency : totalHaplotypeFrequencies )
            {
                frequency /= haplotypeFrequencySum;
            }

            const auto variantQualities =
                VariantQualityCalculator( totalHaplotypeFrequencies, genotypeLikelihoodsAllSamples, readRangesPerSample,
                                          m_samples, variants, mergedHaplotypes, perSampleGenotypes ).getQualities();

            return ModelResults( calledGenotypes, genotypeMetadataPerSample, variantQualities,
                                 variantAnnotationPerSample );
        }

        //-----------------------------------------------------------------------------------------
    }
}
}
