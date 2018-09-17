// All content Copyright (C) 2018 Genomics plc
#ifndef DIPLOIDANNOTATE_HPP
#define DIPLOIDANNOTATE_HPP

#include "utils/matrix.hpp"
#include "stats/functions.hpp"
#include "variant/genotype.hpp"
#include "caller/diploid/genotypeUtils.hpp"
#include "caller/callSet.hpp"
#include "caller/diploid/readSupportAccountant.hpp"
#include "io/readRange.hpp"
#include "caller/metadata.hpp"

namespace echidna
{
namespace caller
{
    namespace annotate
    {
        std::vector< model::VariantMetadata > computeVariantAnnotation(
            const utils::matrix_t & probReadsGivenHaplotypes,
            const io::RegionsReads & readRange,
            const int badReadsWindowSize,
            const std::vector< variant::varPtr_t > & candidateVariants,
            const variant::HaplotypeVector & mergedHaplotypes );

        GenotypeMetadata computeGenotypeMetaData( const std::size_t calledGenotypeIndex,
                                                  const variant::GenotypeVector & genotypes,
                                                  const std::vector< double > & genotypeLikelihoods );

        phred_t computePhaseQuality( const std::size_t calledGenotypeIndex,
                                     const variant::GenotypeVector & genotypes,
                                     const std::vector< double > & genotypeLikelihoods );

        phred_t computeGenotypeQuality( const variant::GenotypeVector & genotypes,
                                        const std::vector< double > & genotypeLikelihoods );

        std::vector< phred_t > get_RR_RA_AA_Likelihoods_as_phred_scores(
            const variant::varPtr_t var,
            const variant::HaplotypeVector & haplotypes,
            const variant::GenotypeVector & genotypes,
            const std::vector< double > & genotypeLikelihoods );

        void annotate( Call & call,
                       const std::size_t varIndex,
                       const std::vector< std::vector< caller::model::VariantMetadata > > & variantAnnotationPerSample,
                       const std::vector< caller::GenotypeMetadata > & genotypeMetadata,
                       const bool outputPhasedGenotypes,
                       const int64_t phaseID );

    }  // namespace annotate
}  // namespace caller
}  // namespace echidna

#endif
