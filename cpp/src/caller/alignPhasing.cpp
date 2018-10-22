// All content Copyright (C) 2018 Genomics plc
#include "caller/alignPhasing.hpp"

namespace wecall
{
namespace caller
{
    std::vector< variant::variantSet_t > generateVariantSetsFromCalls( const callVector_t & calls,
                                                                       const size_t sampleIndex,
                                                                       const size_t ploidy )
    {
        std::vector< variant::variantSet_t > variantsCluster;
        for ( unsigned int ploidyIndex = 0; ploidyIndex < ploidy; ploidyIndex++ )
        {
            std::vector< variant::varPtr_t > variants;
            for ( const auto & call : calls )
            {
                // UNKNOWN state means this variant is not called for this haplotype (treat it as REF)
                if ( call.samples[sampleIndex].genotypeCalls[ploidyIndex] == Call::VAR )
                {
                    variants.emplace_back( call.var );
                }
            }
            variant::variantSet_t variantSet( variants.cbegin(), variants.cend() );
            variantsCluster.emplace_back( variantSet );
        }
        return variantsCluster;
    }

    variant::HaplotypeVector generateHaplotypesFromPhasedClusters(
        const std::vector< variant::variantSet_t > & variantsCluster1,
        const std::vector< variant::variantSet_t > & variantsCluster2,
        const size_t ploidy,
        const utils::referenceSequencePtr_t & referenceSequence,
        const Region & region )
    {
        // generate possible haplotypes spanning both clusters
        variant::HaplotypeVector combinedHaplotypes( region, referenceSequence );
        for ( unsigned int i = 0; i < ploidy; i++ )
        {
            auto variants1 = variantsCluster1[i];
            for ( unsigned int j = 0; j < ploidy; j++ )
            {
                auto variants = variant::variantSet_t( variants1 );
                variants.insert( variantsCluster2[j].cbegin(), variantsCluster2[j].cend() );
                const auto haplotypeId = i * ploidy + j;
                combinedHaplotypes.push_back( variants, haplotypeId );
            }
        }

        // remove duplicate haplotypes as the model cannot handle them
        combinedHaplotypes.sort();
        combinedHaplotypes.merge();
        return combinedHaplotypes;
    }

    bool isHomozygousCluster( const callVector_t calls, const size_t sampleIndex, const size_t ploidy )
    {
        for ( const auto & call : calls )
        {
            const auto expectedValue = call.samples[sampleIndex].genotypeCalls[0];
            for ( unsigned int ploidyIndex = 0; ploidyIndex < ploidy; ploidyIndex++ )
            {
                if ( call.samples[sampleIndex].genotypeCalls[ploidyIndex] != expectedValue )
                {
                    return false;
                }
            }
        }
        return true;
    }

    void alignPhasingOfCallsBasedOnGenotypeAndMergedHaplotypes( const variant::genotypePtr_t & genotype,
                                                                const size_t ploidy,
                                                                const size_t sampleIndex,
                                                                const int64_t phaseIdPrevCluster,
                                                                const variant::HaplotypeVector & combinedHaplotypes,
                                                                const bool cluster2Hom,
                                                                callVector_t & callsCurrentCluster )
    {
        // only align phase for consistent genotype
        // the haplotypeId allows us to backcalculate the index in cluster 1 and cluster2
        std::vector< unsigned int > hap2IndexOrder;
        std::vector< unsigned int > hap1IndexOrder;
        for ( const auto & hapIndex : genotype->getHaplotypeIndices() )
        {
            const size_t haplotypeId = combinedHaplotypes[hapIndex].getId();
            hap1IndexOrder.push_back( haplotypeId / ploidy );
            hap2IndexOrder.push_back( haplotypeId % ploidy );
        }

        // figure out when called genotype is inconsistent
        // - consistent when cluster 2 is homozygous
        // - consistent when cluster 1 and cluster 2 is homozygous (already dealt with earlier)
        // - inconsistent when cluster 1 is homozygous and cluster 2 is heteorozygous
        // - otherwise inconsistent if non-unique haplotypes
        const auto itHap1 = std::unique( hap1IndexOrder.begin(), hap1IndexOrder.end() );
        const auto itHap2 = std::unique( hap2IndexOrder.begin(), hap2IndexOrder.end() );
        const bool unique = ( ( itHap1 == hap1IndexOrder.end() ) and ( itHap2 == hap2IndexOrder.end() ) );

        // try not to align for inconsistent genotype
        if ( not( cluster2Hom or unique ) )
        {
            return;
        }

        // apply phasing and phase set id (ignoring ref calls)
        for ( auto & call : callsCurrentCluster )
        {
            if ( not call.isRefCall() )
            {
                // apply phaseSetId from previous cluster
                auto & phaseSetIdPointer = call.samples[sampleIndex].getAnnotation( Annotation::PS );
                phaseSetIdPointer = phaseIdPrevCluster;

                // shuffle the phase by the order of the c1_index in the genotype
                std::vector< Call::Type > allignedGenotypeCalls;

                for ( unsigned int i = 0; i < ploidy; i++ )
                {
                    allignedGenotypeCalls.push_back(
                        call.samples[sampleIndex].genotypeCalls[hap2IndexOrder[hap1IndexOrder[i]]] );
                }

                for ( unsigned int i = 0; i < ploidy; i++ )
                {
                    call.samples[sampleIndex].genotypeCalls[i] = allignedGenotypeCalls[i];
                }
            }
        }
    }

    void alignPhasingForSample( callVector_t & callsCurrentCluster,
                                const callVector_t & callsPrevCluster,
                                const size_t ploidy,
                                const size_t sampleIndex,
                                const std::string & sampleName,
                                const utils::referenceSequencePtr_t & referenceSequence,
                                const io::perSampleRegionsReads_t & mergedRegionReads,
                                const Region & regionPrevCluster,
                                const Region & combinedRegion )
    {
        // Idea:
        // (1) generate possible combinations of haplotypes from merged calls of the two clusters
        // (2) run model and determine preferred genotype
        // (3) re-annotate the phasing depending on preferred genotype

        const bool cluster1Hom = isHomozygousCluster( callsPrevCluster, sampleIndex, ploidy );
        const bool cluster2Hom = isHomozygousCluster( callsCurrentCluster, sampleIndex, ploidy );
        // do not try to align if first cluster is hom and second het
        if ( cluster1Hom and not cluster2Hom )
        {
            return;
        }

        // get all called variants for sample
        std::vector< variant::varPtr_t > variants;
        for ( auto & call : callsPrevCluster )
        {
            if ( call.samples[sampleIndex].hasVar() )
            {
                variants.push_back( call.var );
            }
        }
        for ( auto & call : callsCurrentCluster )
        {
            if ( call.samples[sampleIndex].hasVar() )
            {
                variants.push_back( call.var );
            }
        }

        // get all variantSets (haplotypes) for both clusters
        std::vector< variant::variantSet_t > variantsCluster1 =
            generateVariantSetsFromCalls( callsPrevCluster, sampleIndex, ploidy );
        std::vector< variant::variantSet_t > variantsCluster2 =
            generateVariantSetsFromCalls( callsCurrentCluster, sampleIndex, ploidy );
        variant::HaplotypeVector combinedHaplotypes = generateHaplotypesFromPhasedClusters(
            variantsCluster1, variantsCluster2, ploidy, referenceSequence, combinedRegion );

        // run model for combined haplotypes for one sample
        model::Model m_model( 1, ploidy * ploidy, {sampleName} );
        const auto results = m_model.getResults( mergedRegionReads, combinedHaplotypes, variants, {ploidy} );

        // change annotation in correspondence to called genotype
        const auto genotypes = results.getCalledGenotypes();

        // only try to align genotype which is well defined
        // a genotype can be nullptr in case of missing read data or genotype size is 0
        if ( genotypes.size() == 1 and genotypes[0] != nullptr )
        {
            const auto genotype = genotypes[0];

            // determine id of previous phase set
            int64_t phaseIdPrevCluster = -1;
            for ( const auto & call : callsPrevCluster )
            {
                if ( not call.isRefCall() )
                {
                    phaseIdPrevCluster = call.samples[sampleIndex].getAnnotation( Annotation::PS );
                    break;
                }
            }

            if ( phaseIdPrevCluster < 0 )
                return;

            // align phases of clusters
            alignPhasingOfCallsBasedOnGenotypeAndMergedHaplotypes( genotype, ploidy, sampleIndex, phaseIdPrevCluster,
                                                                   combinedHaplotypes, cluster2Hom,
                                                                   callsCurrentCluster );
        }
        else
        {
            ECHIDNA_LOG( SUPER_DEBUG, "No genotype data or unexpected number of genotypes: num=" +
                                          std::to_string( genotypes.size() ) + " " + combinedRegion.toString() );
        }
    }

    void alignPhasingBetweenClusters( callVector_t & callsCurrentCluster,
                                      const callVector_t & callsPrevCluster,
                                      const Region & regionCluster,
                                      const Region & regionPrevCluster,
                                      const Region & combinedRegion,
                                      const utils::referenceSequencePtr_t & referenceSequence,
                                      const std::vector< size_t > & ploidyPerSample,
                                      const io::perSampleRegionsReads_t mergedRegionReads,
                                      const std::vector< std::string > & sampleNames )
    {
        // do not try to align if one of the clusters does not contain variants
        if ( callsCurrentCluster.empty() or callsPrevCluster.empty() )
        {
            return;
        }

        // do not try to align if no reads overlap both clusters
        bool readsOverlap = false;
        if ( regionPrevCluster.end() - 1 >= regionCluster.start() + 1 )
        {
            readsOverlap = true;
        }
        else
        {
            // use region of last variant from prev cluster and first variant from current cluster
            // ignore ref calls
            variant::varPtr_t firstVariant = nullptr;
            for ( const auto & call : callsCurrentCluster )
            {
                if ( not call.isRefCall() )
                {
                    firstVariant = call.var;
                    break;
                }
            }

            variant::varPtr_t lastVariant = nullptr;
            for ( const auto & call : boost::adaptors::reverse( callsPrevCluster ) )
            {
                if ( not call.isRefCall() )
                {
                    lastVariant = call.var;
                    break;
                }
            }

            if ( lastVariant == nullptr or firstVariant == nullptr )
            {
                readsOverlap = false;
            }
            else
            {
                readsOverlap =
                    io::readsOverlappingRegions( mergedRegionReads, lastVariant->region(), firstVariant->region() );
            }
        }

        if ( not readsOverlap )
        {
            return;
        }

        // align phase of clusters
        for ( std::size_t sampleIndex = 0; sampleIndex < ploidyPerSample.size(); ++sampleIndex )
        {
            alignPhasingForSample( callsCurrentCluster, callsPrevCluster, ploidyPerSample[sampleIndex], sampleIndex,
                                   sampleNames[sampleIndex], referenceSequence, mergedRegionReads, regionPrevCluster,
                                   combinedRegion );
        }
    }
}
}