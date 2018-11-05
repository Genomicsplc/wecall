// All content Copyright (C) 2018 Genomics plc
#include <iterator>
#include <boost/range/adaptor/reversed.hpp>
#include "caller/diploid/diploid.hpp"
#include "caller/callSet.hpp"
#include "io/readDataSet.hpp"
#include "io/readRange.hpp"
#include "variant/haplotype.hpp"
#include "variant/type/variant.hpp"

#ifndef WECALL_ALIGNPHASING_H
#define WECALL_ALIGNPHASING_H

namespace wecall
{
namespace caller
{

    std::vector< variant::variantSet_t > generateVariantSetsFromCalls( const callVector_t & calls,
                                                                       const size_t sampleIndex,
                                                                       const size_t ploidy );

    variant::HaplotypeVector generateHaplotypesFromPhasedClusters(
        const std::vector< variant::variantSet_t > & variantsCluster1,
        const std::vector< variant::variantSet_t > & variantsCluster2,
        const size_t ploidy,
        const utils::referenceSequencePtr_t & referenceSequence,
        const Region & region );

    bool isHomozygousCluster( const callVector_t calls, const size_t sampleIndex, const size_t ploidy );

    void alignPhasingOfCallsBasedOnGenotypeAndMergedHaplotypes( const variant::genotypePtr_t & genotype,
                                                                const size_t ploidy,
                                                                const size_t sampleIndex,
                                                                const int64_t phaseIdPrevCluster,
                                                                const variant::HaplotypeVector & combinedHaplotypes,
                                                                const bool cluster2Hom,
                                                                callVector_t & callsCurrentCluster );

    void alignPhasingForSample( callVector_t & callsCurrentCluster,
                                const callVector_t & callsPrevCluster,
                                const size_t ploidy,
                                const size_t sampleIndex,
                                const std::string & sampleName,
                                const utils::referenceSequencePtr_t & referenceSequence,
                                const io::perSampleRegionsReads_t & mergedRegionReads,
                                const Region & regionPrevCluster,
                                const Region & region );

    void alignPhasingBetweenClusters( callVector_t & callsCurrentCluster,
                                      const callVector_t & callsPrevCluster,
                                      const Region & regionCluster,
                                      const Region & regionPrevCluster,
                                      const Region & combinedRegion,
                                      const utils::referenceSequencePtr_t & referenceSequence,
                                      const std::vector< size_t > & ploidyPerSample,
                                      const io::perSampleRegionsReads_t mergedRegionReads,
                                      const std::vector< std::string > & sampleNames );
}
}

#endif  // WECALL_ALIGNPHASING_H
