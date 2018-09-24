// All content Copyright (C) 2018 Genomics plc
#ifndef DIPLOID_HPP
#define DIPLOID_HPP

#include "common.hpp"

#include "variant/type/variant.hpp"
#include "variant/haplotype.hpp"
#include "variant/genotype.hpp"
#include "caller/callSet.hpp"
#include "caller/diploid/readSupportAccountant.hpp"
#include "io/read.hpp"
#include "io/readRange.hpp"
#include "io/readDataSet.hpp"
#include "utils/matrix.hpp"
#include "caller/metadata.hpp"

namespace echidna
{
namespace caller
{

    namespace model
    {
        class VCFCallVectorBuilder
        {
        public:
            VCFCallVectorBuilder( const bool makeRefCalls,
                                  const bool outputPhasedCalls,
                                  const bool outputAllVariants,
                                  const bool genotypingMode,
                                  const std::size_t nSamples );

            callVector_t getAnnotatedVariantCalls(
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
                const double referenceCallQualityDeltaThreshold ) const;

        private:
            const bool m_makeRefCalls;
            const bool m_outputPhasedCalls;
            const bool m_outputAllVariants;
            const bool m_genotypingMode;
            const std::size_t m_nSamples;
        };

        class ModelResults
        {
        public:
            ModelResults( const std::vector< variant::genotypePtr_t > & calledGenotypes,
                          const std::vector< GenotypeMetadata > & genotypeMetadata,
                          const std::vector< double > & variantQualities,
                          const std::vector< std::vector< VariantMetadata > > & variantMetadata )
                : m_calledGenotypes( calledGenotypes ),
                  m_genotypeMetadata( genotypeMetadata ),
                  m_variantQualities( variantQualities ),
                  m_variantMetadata( variantMetadata )
            {
            }

            const std::vector< variant::genotypePtr_t > & getCalledGenotypes() const { return m_calledGenotypes; };
            const std::vector< GenotypeMetadata > & getGenotypeMetadata() const { return m_genotypeMetadata; };

            const std::vector< double > & getVariantQualities() const { return m_variantQualities; };
            const std::vector< std::vector< VariantMetadata > > & getVariantMetadata() const
            {
                return m_variantMetadata;
            };

        private:
            const std::vector< variant::genotypePtr_t > m_calledGenotypes;
            const std::vector< GenotypeMetadata > m_genotypeMetadata;

            const std::vector< double > m_variantQualities;
            const std::vector< std::vector< VariantMetadata > > m_variantMetadata;
        };

        struct VariantQualities
        {
            double sumLogProbTotal = constants::unknownValue;
            double sumLogProbNoVariant = constants::unknownValue;
        };

        class Model
        {
        public:
            Model( const int badReadsWindowSize,
                   const std::size_t maxHaplotypesPerCluster,
                   const std::vector< std::string > & samples );

            ModelResults getResults( const io::perSampleRegionsReads_t & readRangesPerSample,
                                     const variant::HaplotypeVector & mergedHaplotypes,
                                     const std::vector< variant::varPtr_t > & variants,
                                     const std::vector< std::size_t > & perSamplePloidy ) const;

        private:
            void computeResultsPerSample( const std::string & sample,
                                          const std::size_t ploidy,
                                          const io::RegionsReads & readRange,
                                          const std::vector< variant::varPtr_t > & variants,
                                          const variant::HaplotypeVector & mergedHaplotypes,
                                          std::map< std::string, variant::GenotypeVector > & perSampleGenotypes,
                                          variant::genotypePtr_t & calledGenotype,
                                          GenotypeMetadata & genotypeMetadata,
                                          std::vector< VariantMetadata > & variantAnnotation,
                                          double & variantQualitiesTotal,
                                          std::vector< double > & reweightedVariantQualities ) const;

            const int m_badReadsWindowSize;
            const std::size_t m_maxHaplotypesPerCluster;
            const std::vector< std::string > m_samples;
        };
    }
}
}

#endif
