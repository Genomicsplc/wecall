// All content Copyright (C) 2018 Genomics plc
#include <algorithm>
#include "caller/metadata.hpp"

#include "caller/diploid/diploidAnnotate.hpp"
#include "io/readUtils.hpp"
#include "io/readRange.hpp"
#include "utils/median.hpp"
#include "stats/models.hpp"

namespace echidna
{
namespace caller
{
    namespace annotate
    {
        std::pair< utils::matrix_t, utils::matrix_t > computePosteriorForReadSupportingVariantsAndReference(
            const utils::matrix_t & probReadsGivenHaplotypes,
            const std::vector< variant::varPtr_t > & candidateVariants,
            const variant::HaplotypeVector & mergedHaplotypes )
        {
            const std::size_t nReads = probReadsGivenHaplotypes.size1();

            // compute denominator
            std::vector< double > sumOfReadLikelihoodsOverAllHaplotypes( nReads );
            for ( std::size_t readIndex = 0; readIndex < nReads; ++readIndex )
            {
                const utils::matrixRow_t probReadGivenHaplotypes( probReadsGivenHaplotypes, readIndex );
                sumOfReadLikelihoodsOverAllHaplotypes[readIndex] =
                    utils::sumMatrixRowOverAllIndices( probReadGivenHaplotypes );
            }

            // compute matrix p(v|r) and p(ref_at_v|r)
            const auto nVariants = candidateVariants.size();
            utils::matrix_t resultsForVariants( nVariants, nReads );
            utils::matrix_t resultsForReference( nVariants, nReads );

            for ( std::size_t variantIndex = 0; variantIndex < nVariants; ++variantIndex )
            {
                const auto thisVar = candidateVariants[variantIndex];

                const auto & haplotypeIndicesForThisVariant = mergedHaplotypes.getHaplotypeIndicesForVariant( thisVar );

                const auto & haplotypeIndicesForReferenceAtThisVariantPosition =
                    mergedHaplotypes.getIndicesForReference( thisVar->region() );

                for ( std::size_t readIndex = 0; readIndex < nReads; ++readIndex )
                {
                    const utils::matrixRow_t probReadGivenHaplotypes( probReadsGivenHaplotypes, readIndex );

                    resultsForVariants( variantIndex, readIndex ) =
                        utils::sumMatrixRowOverIndexSubset( probReadGivenHaplotypes, haplotypeIndicesForThisVariant );

                    resultsForReference( variantIndex, readIndex ) = utils::sumMatrixRowOverIndexSubset(
                        probReadGivenHaplotypes, haplotypeIndicesForReferenceAtThisVariantPosition );
                }
            }

            // normalise
            for ( std::size_t readIndex = 0; readIndex < nReads; ++readIndex )
            {

                const auto denominator = sumOfReadLikelihoodsOverAllHaplotypes[readIndex];
                if ( denominator > std::numeric_limits< double >::min() )
                {
                    for ( std::size_t variantIndex = 0; variantIndex < nVariants; ++variantIndex )
                    {
                        resultsForVariants( variantIndex, readIndex ) /= denominator;
                        resultsForReference( variantIndex, readIndex ) /= denominator;
                    }
                }
                else
                {
                    for ( std::size_t variantIndex = 0; variantIndex < nVariants; ++variantIndex )
                    {
                        resultsForVariants( variantIndex, readIndex ) = 0.0;
                        resultsForReference( variantIndex, readIndex ) = 0.0;
                    }
                }
            }

            return std::make_pair( resultsForVariants, resultsForReference );
        }

        std::vector< model::VariantMetadata > computeVariantAnnotation(
            const utils::matrix_t & probReadsGivenHaplotypes,
            const io::RegionsReads & readRange,
            const int badReadsWindowSize,
            const std::vector< variant::varPtr_t > & candidateVariants,
            const variant::HaplotypeVector & mergedHaplotypes )
        {
            const auto nVariants = candidateVariants.size();

            std::vector< model::VariantMetadata > varAnnotations;
            for ( std::size_t variantIndex = 0; variantIndex < nVariants; ++variantIndex )
            {
                const auto variantRegions = candidateVariants[variantIndex]->getStartEndRegions(
                    mergedHaplotypes.region().start(), mergedHaplotypes.region().end() );

                varAnnotations.emplace_back( variantRegions, badReadsWindowSize );
            }

            const auto posteriorForReadSupportingVariantsAndReference =
                computePosteriorForReadSupportingVariantsAndReference( probReadsGivenHaplotypes, candidateVariants,
                                                                       mergedHaplotypes );
            const auto & posteriorForReadSupportingVariants = posteriorForReadSupportingVariantsAndReference.first;
            const auto & posteriorForReadSupportingReference = posteriorForReadSupportingVariantsAndReference.second;

            for ( std::size_t variantIndex = 0; variantIndex < nVariants; ++variantIndex )
            {
                std::size_t readIndex = 0;
                const auto variantRegions = varAnnotations[variantIndex].regions();

                for ( const auto & read : readRange )
                {
                    const auto readInterval = read.getMaximalReadInterval();

                    bool considerRead = false;
                    for ( const auto & region : variantRegions )
                    {
                        if ( readInterval.contains( region.interval() ) and readInterval.overlaps( region.interval() ) )
                        {
                            considerRead = true;
                        }
                    }

                    // This relies too much on the original position in the bam. We should be able to infer this.
                    if ( considerRead )
                    {
                        const double probReadSupportsThisVar =
                            posteriorForReadSupportingVariants( variantIndex, readIndex );

                        const double probReadSupportsReference =
                            posteriorForReadSupportingReference( variantIndex, readIndex );

                        varAnnotations[variantIndex].accountForRead( probReadSupportsThisVar, probReadSupportsReference,
                                                                     read );
                    }

                    ++readIndex;
                }
            }
            return varAnnotations;
        }

        GenotypeMetadata computeGenotypeMetaData( const std::size_t calledGenotypeIndex,
                                                  const variant::GenotypeVector & genotypes,
                                                  const std::vector< double > & genotypeLikelihoods )
        {
            GenotypeMetadata genotypeMetadata;

            genotypeMetadata.phaseQuality =
                annotate::computePhaseQuality( calledGenotypeIndex, genotypes, genotypeLikelihoods );

            genotypeMetadata.genotypeQuality = annotate::computeGenotypeQuality( genotypes, genotypeLikelihoods );

            return genotypeMetadata;
        }

        phred_t computePhaseQuality( const std::size_t calledGenotypeIndex,
                                     const variant::GenotypeVector & genotypes,
                                     const std::vector< double > & genotypeLikelihoods )
        {
            assert( genotypes.size() == genotypeLikelihoods.size() );

            double maxLikelihood = 0.0;
            double errorLikelihood = 0.0;

            for ( const auto & genotypeIndex :
                  genotypes.genotypesWithSameNonPhasedRepresentation( calledGenotypeIndex ) )
            {
                const double likelihood =
                    genotypeLikelihoods[genotypeIndex] * genotypes[genotypeIndex]->nCombinationsThisGenotype();

                if ( likelihood > maxLikelihood )
                {
                    errorLikelihood += maxLikelihood;
                    maxLikelihood = likelihood;
                }
                else
                {
                    errorLikelihood += likelihood;
                }
            }

            // Return probability of the called genotype as a Phred score
            return stats::toPhredQ( errorLikelihood / ( maxLikelihood + errorLikelihood ) );
        }

        phred_t computeGenotypeQuality( const variant::GenotypeVector & genotypes,
                                        const std::vector< double > & genotypeLikelihoods )
        {
            double maxLikelihood = 0.0;
            double errorLikelihood = 0.0;

            for ( std::size_t genotypeIndex = 0; genotypeIndex < genotypes.size(); ++genotypeIndex )
            {
                const double likelihood =
                    genotypeLikelihoods[genotypeIndex] * genotypes[genotypeIndex]->nCombinationsThisGenotype();

                if ( likelihood > maxLikelihood )
                {
                    errorLikelihood += maxLikelihood;
                    maxLikelihood = likelihood;
                }
                else
                {
                    errorLikelihood += likelihood;
                }
            }

            // Return probability of the called genotype as a Phred score
            return stats::toPhredQ( errorLikelihood / ( maxLikelihood + errorLikelihood ) );
        }

        std::vector< phred_t > get_RR_RA_AA_Likelihoods_as_phred_scores(
            const variant::varPtr_t var,
            const variant::HaplotypeVector & haplotypes,
            const variant::GenotypeVector & genotypes,
            const std::vector< double > & genotypeLikelihoods )
        {
            std::vector< double > totalLikelihood( 1 + genotypes.ploidy(), 0.0 );  // RR, RA, AA
            std::vector< phred_t > retVals( 1 + genotypes.ploidy() );

            for ( std::size_t genotypeIndex = 0; genotypeIndex < genotypes.size(); ++genotypeIndex )
            {
                const auto nStrandsWithVar = genotypes[genotypeIndex]->nStrandsContainingVariant( haplotypes, var );
                const auto factor = genotypes[genotypeIndex]->nCombinationsThisGenotype();

                totalLikelihood[nStrandsWithVar] += genotypeLikelihoods[genotypeIndex] * factor;
            }

            std::transform( totalLikelihood.begin(), totalLikelihood.end(), totalLikelihood.begin(), log10 );

            caller::model::rescaleLogLikelihoods( totalLikelihood, true );

            const auto toPhred = []( double val ) -> phred_t
            {
                return std::min( stats::phredCoefficient * val, constants::maxPhredScore );
            };

            std::transform( totalLikelihood.begin(), totalLikelihood.end(), retVals.begin(), toPhred );

            return retVals;
        }

        template < typename T >
        T bespokeMaxElement( const std::vector< T > & values )
        {
            const auto elem = std::max_element( values.begin(), values.end() );
            return elem == values.end() ? std::numeric_limits< T >::quiet_NaN() : *elem;
        }

        void annotate( Call & call,
                       const std::size_t varIndex,
                       const std::vector< std::vector< caller::model::VariantMetadata > > & variantAnnotationPerSample,
                       const std::vector< caller::GenotypeMetadata > & genotypeMetadata,
                       const bool outputPhasedGenotypes,
                       const int64_t phaseID )
        {
            // Add quality as an annotation, so we can see it for each variant when those at
            // the same position are combined.

            call.addAnnotation( Annotation::PP, call.qual );

            // Initialise coverage totals and their "Where Called" (WC) equivalents
            int64_t nTotalReadDepthAllSamples = 0, nTotalReverseReadDepthAllSamples = 0,
                    nTotalForwardReadDepthAllSamples = 0;
            int64_t TC_WC = 0, TCR_WC = 0, TCF_WC = 0, VC_WC = 0, VCR_WC = 0, VCF_WC = 0;

            for ( std::size_t sampleIndex = 0; sampleIndex < call.samples.size(); ++sampleIndex )
            {
                call.samples[sampleIndex].addAnnotation(
                    Annotation::PL, variantAnnotationPerSample[sampleIndex][varIndex].genotypeLikelihoods );
            }

            for ( std::size_t sampleIndex = 0; sampleIndex < call.samples.size(); ++sampleIndex )
            {
                call.samples[sampleIndex].addAnnotation( Annotation::GQ,
                                                         genotypeMetadata[sampleIndex].genotypeQuality );
            }

            if ( outputPhasedGenotypes )
            {
                for ( std::size_t sampleIndex = 0; sampleIndex < call.samples.size(); ++sampleIndex )
                {
                    call.samples[sampleIndex].addAnnotation( Annotation::PQ,
                                                             genotypeMetadata[sampleIndex].phaseQuality );

                    call.samples[sampleIndex].addAnnotation( Annotation::PS, phaseID );
                }
            }

            // Genotype (sample) specific annotations.
            std::vector< double > rootMeanSquareMappingQuals = {};
            std::vector< phred_t > medianBaseQuals = {};
            for ( std::size_t sampleIndex = 0; sampleIndex < call.samples.size(); ++sampleIndex )
            {
                if ( call.samples[sampleIndex].hasVar() )
                {
                    const auto medianMinBaseQualities =
                        variantAnnotationPerSample[sampleIndex][varIndex].getMedianMinBaseQualities();

                    if ( not std::isnan( medianMinBaseQualities ) )
                    {
                        medianBaseQuals.push_back( medianMinBaseQualities );
                    }

                    const auto rootMeanSquareMappingQuality =
                        variantAnnotationPerSample[sampleIndex][varIndex].getRootMeanSquareMappingQuality();
                    if ( not std::isnan( rootMeanSquareMappingQuality ) )
                    {
                        rootMeanSquareMappingQuals.push_back( rootMeanSquareMappingQuality );
                    }
                }
            }

            // call sample annotations
            // Contains all per-(variant & sample) information gained from reading the input file.
            int64_t nForwardReadsSupportingVariant = 0;
            int64_t nReverseReadsSupportingVariant = 0;
            for ( std::size_t sampleIndex = 0; sampleIndex < call.samples.size(); ++sampleIndex )
            {
                const auto sampleReadAccountant = variantAnnotationPerSample[sampleIndex][varIndex];

                const auto nVariantSupportingReverseReads = sampleReadAccountant.getVariantSupportingReverseReads();
                const auto nVariantSupportingForwardReads = sampleReadAccountant.getVariantSupportingForwardReads();
                const auto nVariantSupportingReads = nVariantSupportingForwardReads + nVariantSupportingReverseReads;

                const auto nRefSupportingReads = sampleReadAccountant.getReferenceSupportingReads();

                const auto nTotalReverseReadDepth = sampleReadAccountant.getTotalReverseReads();
                const auto nTotalForwardReadDepth = sampleReadAccountant.getTotalForwardReads();
                const auto nTotalReadDepth = nTotalReverseReadDepth + nTotalForwardReadDepth;

                const bool isCalledInThisSample = call.samples[sampleIndex].hasVar();

                nTotalReadDepthAllSamples += nTotalReadDepth;
                nTotalReverseReadDepthAllSamples += nTotalReverseReadDepth;
                nTotalForwardReadDepthAllSamples += nTotalForwardReadDepth;

                if ( isCalledInThisSample )
                {
                    TC_WC += nTotalReadDepth;
                    TCR_WC += nTotalReverseReadDepth;
                    TCF_WC += nTotalForwardReadDepth;

                    VC_WC += nVariantSupportingReads;
                    VCR_WC += nVariantSupportingReverseReads;
                    VCF_WC += nVariantSupportingForwardReads;
                }

                nForwardReadsSupportingVariant += nVariantSupportingForwardReads;
                nReverseReadsSupportingVariant += nVariantSupportingReverseReads;

                const std::vector< int64_t > allelicDepths = {nRefSupportingReads, nVariantSupportingReads};
                call.samples[sampleIndex].addAnnotation( Annotation::AD, allelicDepths );
                call.samples[sampleIndex].addAnnotation( Annotation::FORMAT_DP, nTotalReadDepth );

                const double allelicFrequency =
                    nTotalReadDepth == 0
                        ? std::numeric_limits< double >::quiet_NaN()
                        : static_cast< double >( nVariantSupportingReads ) / static_cast< double >( nTotalReadDepth );

                call.samples[sampleIndex].addAnnotation( Annotation::VAF, {allelicFrequency} );
            }

            const auto nReadsSupportingVariant = nForwardReadsSupportingVariant + nReverseReadsSupportingVariant;

            call.addAnnotation( Annotation::DP, nTotalReadDepthAllSamples );
            call.addAnnotation( Annotation::DPR, nTotalReverseReadDepthAllSamples );
            call.addAnnotation( Annotation::DPF, nTotalForwardReadDepthAllSamples );

            call.addAnnotation( Annotation::VC, nReadsSupportingVariant );
            call.addAnnotation( Annotation::VCR, nReverseReadsSupportingVariant );
            call.addAnnotation( Annotation::VCF, nForwardReadsSupportingVariant );

            call.addAnnotation( Annotation::ABPV, stats::probSufficientVarCoverageToSupportHet( TC_WC, VC_WC ) );
            call.addAnnotation( Annotation::SBPV,
                                stats::probVarSupportNotBiasedByStrand( TCF_WC, VCF_WC, TCR_WC, VCR_WC ) );

            call.addAnnotation( Annotation::MQ, bespokeMaxElement( rootMeanSquareMappingQuals ) );

            call.addAnnotation( Annotation::BR, bespokeMaxElement( medianBaseQuals ) );
            call.addAnnotation( Annotation::QD, stats::variantSupportPerRead( call.var->phredScaledPrior(), call.qual,
                                                                              nReadsSupportingVariant ) );
        }

    }  // namespace annotate
}  // namespace caller
}  // namespace echidna
