// All content Copyright (C) 2018 Genomics plc
#ifndef BREAKPOINT_VARIANT_GENERATOR
#define BREAKPOINT_VARIANT_GENERATOR

#include <io/readRange.hpp>
#include "utils/referenceSequence.hpp"
#include "assembly/sequenceGraph.hpp"
#include "variant/type/breakpoint.hpp"
#include "variant/type/variant.hpp"

namespace echidna
{
namespace variant
{
    class BreakpointVariantGenerator
    {
    public:
        BreakpointVariantGenerator( const utils::referenceSequencePtr_t & referenceSequence,
                                    const std::size_t kmerSize,
                                    const std::size_t maxKmerSize,
                                    const std::size_t minChainSupport = 1 );

        variantSet_t getVariantCandidates( const breakpointLocusSet_t & breakpoints,
                                           const io::perSampleRegionsReads_t & reads ) const;

    private:
        caller::Region referenceRegionFromBreakpoint( const breakpointLocusPtr_t & bp ) const;

        assembly::SequenceGraph processData( const caller::SetRegions & referenceRegions,
                                             const io::perSampleRegionsReads_t & reads ) const;

        std::pair< assembly::SequenceGraph, bool > createSequenceGraph( const caller::SetRegions & referenceRegions,
                                                                        const io::perSampleRegionsReads_t & reads,
                                                                        const std::size_t kmerSize,
                                                                        const bool disallowRepeats ) const;

        const std::size_t m_defaultKmerSize;
        const std::size_t m_maxKmerSize;
        const std::size_t m_minChainSupport;
        const utils::referenceSequencePtr_t m_refSeq;
    };

    class BreakpointClusterer
    {
    public:
        BreakpointClusterer( const std::size_t minDistance ) : m_minDistance( minDistance ) {}

        std::vector< std::vector< breakpointLocusPtr_t > > getClusters(
            const breakpointLocusSet_t & breakpointSet ) const;

    private:
        const std::size_t m_minDistance;
    };

    class MultiLocusBreakpointVariantGenerator
    {
    public:
        MultiLocusBreakpointVariantGenerator( const caller::Region & region,
                                              const utils::referenceSequencePtr_t & referenceSequence,
                                              const std::size_t defaultKmerSize,
                                              const std::size_t maxKmerSize,
                                              const int largeVariantSizeDefinition );

        variantSet_t getVariantCandidates( const breakpointLocusSet_t & breakpointSet,
                                           const io::perSampleRegionsReads_t & allReads ) const;

    private:
        const caller::Region m_region;
        const utils::referenceSequencePtr_t m_referenceSequence;
        const std::size_t m_defaultKmerSize;
        const std::size_t m_maxKmerSize;
        const int m_largeVariantSizeDefinition;
    };
}
}

#endif
