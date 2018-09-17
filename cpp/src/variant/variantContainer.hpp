// All content Copyright (C) 2018 Genomics plc
#ifndef VARIANT_CONTAINER_HPP
#define VARIANT_CONTAINER_HPP

#include "variant/type/variant.hpp"
#include "variant/type/breakpoint.hpp"

#include <memory>
#include <map>
#include <io/readRange.hpp>

namespace echidna
{
namespace io
{
    class ReadDataset;
    using readDataset_t = std::shared_ptr< const ReadDataset >;
}

namespace variant
{

    /// VariantContainer Class that contains variants along with initial data computed from the bamfiles.
    /// These are stored per variant in nested class VariantCounts
    class VariantContainer
    {

    public:
        /// VariantCounts Nested class that stores per variant & sample information
        // about variant. These are computed from pre-aligned reads.
        struct VariantCounts
        {
            int64_t getPercentVariantCoverage() const;

            int64_t m_totalReads = 0;
            int64_t m_totalVariantSupportingReads = 0;
        };

    public:
        VariantContainer( phred_t minBaseQual, phred_t minMappingQual )
            : m_minBaseQual( minBaseQual ), m_minMappingQual( minMappingQual )
        {
        }

        void addCandidateVariant( const varPtr_t & candidate, const double prior );
        void addGenotypingVariant( const varPtr_t & variant );

        variantSet_t getVariants() const;
        breakpointLocusSet_t getBreakpoints() const;

        void addVariantsFromRead( io::readPtr_t readPtr,
                                  const std::vector< varPtr_t > & variants,
                                  const std::vector< breakpointPtr_t > & breakpoints,
                                  const std::string & sampleName );

        // TODO: Should this be in this class? Coverage for a region is not a property of a variant container.
        void computeCoverage( const caller::Region & blockRegion, io::perSampleRegionsReads_t allReads );

        int64_t totalReadsSupportingVariant( varPtr_t varPtr ) const;
        int64_t maxReadPercentVariantCoverage( varPtr_t varPtr ) const;

    private:
        void addBreakpoints( const std::vector< breakpointPtr_t > & breakpoints );
        void addVariant( varPtr_t var, phred_t baseQual, io::readPtr_t readPtr, const std::string & sampleName );

        std::map< varPtr_t, std::map< std::string, VariantCounts >, varPtrComp > m_variantCountsPerSample;

        std::map< int64_t, breakpointLocusPtr_t > m_startBreakpointLoci;
        std::map< int64_t, breakpointLocusPtr_t > m_endBreakpointLoci;

        const phred_t m_minBaseQual;
        const phred_t m_minMappingQual;
    };
}
}

#endif
