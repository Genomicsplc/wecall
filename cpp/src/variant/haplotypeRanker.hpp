// All content Copyright (C) 2018 Genomics plc
#ifndef HAPLOTYPE_RANKER_HPP
#define HAPLOTYPE_RANKER_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <io/readRange.hpp>

#include "variant/haplotype.hpp"

namespace wecall
{
namespace io
{
    class Read;
}

namespace variant
{
    class AlignmentHaplotypeRanker
    {
    public:
        AlignmentHaplotypeRanker( const io::perSampleRegionsReads_t & reads ) : m_reads( reads ) {}

        std::set< std::size_t > getTopHaplotypes( const HaplotypeVector & haplotypes,
                                                  const uint64_t maxHaplotypes ) const;

    private:
        io::perSampleRegionsReads_t m_reads;
    };
}
}

#endif
