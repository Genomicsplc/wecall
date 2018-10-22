// All content Copyright (C) 2018 Genomics plc
#ifndef ECHIDNA_VARIANT_QUALITY_CALCULATOR_HPP
#define ECHIDNA_VARIANT_QUALITY_CALCULATOR_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include "utils/logging.hpp"
#include "io/readRange.hpp"
#include "variant/genotype.hpp"

namespace wecall
{
namespace caller
{
    namespace model
    {
        class VariantQualityCalculator
        {
        public:
            VariantQualityCalculator() = delete;

            VariantQualityCalculator( const std::vector< double > & haplotypeFrequencies,
                                      const std::vector< std::vector< double > > & genotypeLikelihoods,
                                      const io::perSampleRegionsReads_t & readRangesPerSample,
                                      const std::vector< std::string > & samples,
                                      const std::vector< variant::varPtr_t > & candidateVariants,
                                      const variant::HaplotypeVector & candiateHaplotypes,
                                      const std::map< std::string, variant::GenotypeVector > & candidateGenotypes )
                : m_haplotypeFrequencies( haplotypeFrequencies ),
                  m_genotypeLikelihoods( genotypeLikelihoods ),
                  m_readRangesPerSample( readRangesPerSample ),
                  m_samples( samples ),
                  m_candidateVariants( candidateVariants ),
                  m_candidateHaplotypes( candiateHaplotypes ),
                  m_candidateGenotypes( candidateGenotypes )
            {
            }

            std::vector< double > getQualities() const;

        private:
            double computeVariantQuality( const std::size_t variantIndex ) const;

            std::vector< double > computeReweightedHaplotypeFrequencies( const std::size_t variantIndex ) const;

        private:
            const std::vector< double > & m_haplotypeFrequencies;
            const std::vector< std::vector< double > > & m_genotypeLikelihoods;
            const io::perSampleRegionsReads_t & m_readRangesPerSample;
            const std::vector< std::string > m_samples;
            const std::vector< variant::varPtr_t > & m_candidateVariants;
            const variant::HaplotypeVector & m_candidateHaplotypes;
            const std::map< std::string, variant::GenotypeVector > & m_candidateGenotypes;
        };
    }
}
}

#endif
