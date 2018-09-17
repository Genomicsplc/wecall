// All content Copyright (C) 2018 Genomics plc
#ifndef READSUPPORTACCOUNTANT_HPP
#define READSUPPORTACCOUNTANT_HPP

#include "common.hpp"

#include "utils/logging.hpp"
#include "stats/models.hpp"
#include "io/read.hpp"
#include "caller/region.hpp"

namespace echidna
{
namespace caller
{
    namespace model
    {
        class VariantMetadata
        {
        public:
            VariantMetadata( const caller::SetRegions & regions, const int badReadsWindowSize )
                : m_regions( regions ), m_badReadsWindowSize( badReadsWindowSize )
            {
            }
            void accountForRead( const double probReadSupportsVariant,
                                 const double probReadSupportsReference,
                                 const io::Read & read );

            std::string toString() const;

            int64_t getTotalForwardReads() const;
            int64_t getTotalReverseReads() const;
            int64_t getVariantSupportingForwardReads() const;
            int64_t getVariantSupportingReverseReads() const;
            int64_t getReferenceSupportingForwardReads() const;
            int64_t getReferenceSupportingReverseReads() const;
            int64_t getReferenceSupportingReads() const;

            double getRootMeanSquareMappingQuality() const { return stats::rootMeanSquare( m_mappingQuals ); }
            phred_t getMedianMinBaseQualities() const;

            const caller::SetRegions & regions() const { return m_regions; }

        private:
            const caller::SetRegions m_regions;
            const int m_badReadsWindowSize;

            double m_forwardSupportingVariant = 0.0;
            double m_reverseSupportingVariant = 0.0;

            double m_forwardNotSupportingVariant = 0.0;
            double m_reverseNotSupportingVariant = 0.0;

            double m_forwardSupportingReference = 0.0;
            double m_reverseSupportingReference = 0.0;

            std::vector< phred_t > m_minBaseQualitiesPerRead;
            std::vector< int64_t > m_mappingQuals;

        public:
            std::vector< phred_t > genotypeLikelihoods;
        };
    }
}
}

#endif  // READSUPPORTACCOUNTANT_HPP
