// All content Copyright (C) 2018 Genomics plc
#ifndef GENOTYPE_UTILS_HPP
#define GENOTYPE_UTILS_HPP

#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "utils/matrix.hpp"
#include "variant/genotype.hpp"

namespace echidna
{
namespace caller
{
    namespace model
    {
        template < typename Sequence >
        void rescaleLogLikelihoods( Sequence & sequence, bool hasMultipleReads )
        {
            if ( hasMultipleReads )
            {
                const auto maxLogLikelihood = *std::max_element( sequence.begin(), sequence.end() );
                for ( auto it = sequence.begin(); it != sequence.end(); ++it )
                {
                    *it = log( std::max( std::numeric_limits< double >::min(), exp( *it - maxLogLikelihood ) ) );
                }
            }
            else
            {
                for ( auto it = sequence.begin(); it != sequence.end(); ++it )
                {
                    *it = 1.0;
                }
            }
        }

        std::vector< double > computeGenotypeLikelihoods( const variant::GenotypeVector & genotypes,
                                                          const utils::matrix_t & probReadsGivenHaplotypes,
                                                          const variant::HaplotypeVector & haplotypes );
    }
}
}

#endif
