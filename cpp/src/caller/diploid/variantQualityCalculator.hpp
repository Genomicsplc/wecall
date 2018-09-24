// All content Copyright (C) 2018 Genomics plc
#ifndef ECHIDNA_VARIANT_QUALITY_CALCULATOR_HPP
#define ECHIDNA_VARIANT_QUALITY_CALCULATOR_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include "utils/logging.hpp"
#include "io/readRange.hpp"
#include "variant/genotype.hpp"
#include "diploid.hpp"

namespace echidna
{
namespace caller
{
    namespace model
    {

        class VariantQualityCalculator
        {
        public:
            VariantQualityCalculator() = delete;

            VariantQualityCalculator( const std::vector< std::string > & samples,
                                      const std::vector< variant::varPtr_t > & candidateVariants,
                                      const std::vector< std::vector< double > > & reweightedVariantQualities,
                                      const std::vector< double > & variantQualitiesTotal )
                : m_samples( samples ),
                  m_candidateVariants( candidateVariants ),
                  m_reweightedVariantQualities( reweightedVariantQualities ),
                  m_variantQualitiesTotal( variantQualitiesTotal )
            {
            }

            std::vector< double > getQualities() const;

        private:
            double computeVariantQuality( const std::size_t variantIndex ) const;

        private:
            const std::vector< std::string > m_samples;
            const std::vector< variant::varPtr_t > & m_candidateVariants;
            const std::vector< std::vector< double > > & m_reweightedVariantQualities;
            const std::vector< double > & m_variantQualitiesTotal;
        };
    }
}
}

#endif
