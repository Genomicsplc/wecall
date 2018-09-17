// All content Copyright (C) 2018 Genomics plc
#include <algorithm>
#include <caller/haplotypeLikelihoods.hpp>
#include <boost/algorithm/string/join.hpp>
#include "utils/bestScoreSelector.hpp"
#include "variant/haplotypeRanker.hpp"
#include "variant/haplotype.hpp"
#include "utils/combinations.hpp"
#include "mapping/hashMapper.hpp"
#include "io/read.hpp"

namespace echidna
{
namespace variant
{
    std::set< std::size_t > AlignmentHaplotypeRanker::getTopHaplotypes( const HaplotypeVector & haplotypes,
                                                                        const uint64_t maxHaplotypes ) const
    {
        if ( haplotypes.size() <= maxHaplotypes )
        {
            std::vector< std::size_t > v( haplotypes.size() );
            std::iota( v.begin(), v.end(), 0 );
            return std::set< std::size_t >( v.begin(), v.end() );
        }

        std::vector< double > totalHaplotypeFrequencies( haplotypes.size(), 0.0 );
        for ( const auto & readRangePair : m_reads )
        {
            const auto haplotypeLikelihoods = caller::computeHaplotypeLikelihoods( haplotypes, readRangePair.second );
            const auto haplotypeFrequencies = caller::computeHaplotypeFrequencies( haplotypeLikelihoods );

            for ( std::size_t haplotypeIndex = 0; haplotypeIndex != haplotypes.size(); ++haplotypeIndex )
            {
                totalHaplotypeFrequencies[haplotypeIndex] += haplotypeFrequencies[haplotypeIndex];
            }
        }
        const auto bestIndicies =
            utils::indiciesWithHighestValues< double >( totalHaplotypeFrequencies, maxHaplotypes, 1.0, 1.0e5 );
        return std::set< std::size_t >( bestIndicies.cbegin(), bestIndicies.cend() );
    }

    //-------------------------------------------------------------------------------------------------
}
}
