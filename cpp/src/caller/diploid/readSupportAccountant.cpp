// All content Copyright (C) 2018 Genomics plc
#include <sstream>
#include <cmath>

#include "caller/diploid/readSupportAccountant.hpp"

#include "utils/median.hpp"
#include "io/readUtils.hpp"

namespace wecall
{
namespace caller
{
    namespace model
    {
        //-------------------------------------------------------------------------------------

        void VariantMetadata::accountForRead( const double probReadSupportsVariant,
                                              const double probReadSupportsReference,
                                              const io::Read & read )
        {
            if ( read.isReverse() )
            {
                m_reverseSupportingVariant += probReadSupportsVariant;
                m_reverseNotSupportingVariant += 1.0 - probReadSupportsVariant;
                m_reverseSupportingReference += probReadSupportsReference;
            }
            else
            {
                m_forwardSupportingVariant += probReadSupportsVariant;
                m_forwardNotSupportingVariant += 1.0 - probReadSupportsVariant;
                m_forwardSupportingReference += probReadSupportsReference;
            }

            if ( probReadSupportsVariant >= constants::thresholdForReadSupportingVariant )
            {
                // TODO(ES): What interval should be used for large deletion?

                const auto interval = m_regions.getSpan().interval();

                m_minBaseQualitiesPerRead.push_back(
                    io::read::minBaseQualityInReadAroundInterval( read, interval, m_badReadsWindowSize ) );
                m_mappingQuals.push_back( read.getMappingQuality() );
            }
        }

        phred_t VariantMetadata::getMedianMinBaseQualities() const
        {
            if ( m_minBaseQualitiesPerRead.empty() )
            {
                return std::numeric_limits< phred_t >::quiet_NaN();
            }
            else
            {
                return utils::functional::median( m_minBaseQualitiesPerRead );
            }
        }

        int64_t VariantMetadata::getTotalForwardReads() const
        {
            return static_cast< int64_t >( std::round( m_forwardSupportingVariant + m_forwardNotSupportingVariant ) );
        }

        int64_t VariantMetadata::getTotalReverseReads() const
        {
            return static_cast< int64_t >( std::round( m_reverseSupportingVariant + m_reverseNotSupportingVariant ) );
        }

        int64_t VariantMetadata::getVariantSupportingForwardReads() const
        {
            return static_cast< int64_t >( std::round( m_forwardSupportingVariant ) );
        }

        int64_t VariantMetadata::getVariantSupportingReverseReads() const
        {
            return static_cast< int64_t >( std::round( m_reverseSupportingVariant ) );
        }

        int64_t VariantMetadata::getReferenceSupportingForwardReads() const
        {
            return static_cast< int64_t >( std::round( m_forwardSupportingReference ) );
        }

        int64_t VariantMetadata::getReferenceSupportingReverseReads() const
        {
            return static_cast< int64_t >( std::round( m_reverseSupportingReference ) );
        }

        int64_t VariantMetadata::getReferenceSupportingReads() const
        {
            return this->getReferenceSupportingForwardReads() + this->getReferenceSupportingReverseReads();
        }

        //-------------------------------------------------------------------------------------

        std::string VariantMetadata::toString() const
        {
            std::stringstream ss;
            ss << "Forward supporting variant: " << m_forwardSupportingVariant << "\n"
               << "Reverse supporting variant: " << m_reverseSupportingVariant << "\n"
               << "Forward not supporting variant: " << m_forwardNotSupportingVariant << "\n"
               << "Reverse not supporting variant: " << m_reverseNotSupportingVariant;

            return ss.str();
        }
    }
}
}
