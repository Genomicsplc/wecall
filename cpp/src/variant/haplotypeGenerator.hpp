// All content Copyright (C) 2018 Genomics plc
#ifndef HAPLOTYPE_GENERATOR_HPP
#define HAPLOTYPE_GENERATOR_HPP

#include <vector>

#include "variant/haplotypeRanker.hpp"
#include "variant/haplotype.hpp"
#include "io/readRange.hpp"
#include "io/fastaFile.hpp"

namespace wecall
{
namespace variant
{

    class AlignmentHaplotypeGenerator
    {
    public:
        AlignmentHaplotypeGenerator( const std::vector< variant::varPtr_t > & variants,
                                     const caller::SetRegions & region,
                                     const io::perSampleRegionsReads_t & readsPerSample,
                                     utils::referenceSequencePtr_t referenceSequence,
                                     const int64_t maxHaplotypesPerRanker,
                                     const std::size_t minReadsToSupportClaim );

        HaplotypeVector generateHaplotypes() const;

    protected:
        HaplotypeVector generateRawHaplotypes() const;

        std::vector< variantSet_t > bestVariantCombos( const VariantCluster & cluster ) const;

    private:
        const std::vector< variant::varPtr_t > m_vars;
        const caller::SetRegions m_regions;
        const io::perSampleRegionsReads_t & m_readsPerSample;
        const utils::referenceSequencePtr_t m_referenceSequence;
        const int64_t m_maxHaplotypesPerRanker;
        const std::size_t m_minReadsToSupportClaim;
    };
}
}

#endif
