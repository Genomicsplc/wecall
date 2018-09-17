// All content Copyright (C) 2018 Genomics plc
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "caller/diploid/genotypeUtils.hpp"
#include "utils/indexedProduct.hpp"

namespace echidna
{
namespace caller
{
    namespace model
    {
        std::vector< double > computeGenotypeLikelihoods( const variant::GenotypeVector & genotypes,
                                                          const utils::matrix_t & probReadsGivenHaplotypes,
                                                          const variant::HaplotypeVector & haplotypes )
        {
            std::vector< double > genotypeLikelihoods( genotypes.size() );
            const auto nReads = probReadsGivenHaplotypes.size1();
            ECHIDNA_ASSERT( nReads > 0, "Only makes sense to compute likelihoods if there is read-data" );

            for ( std::size_t genotypeIndex = 0; genotypeIndex < genotypes.size(); ++genotypeIndex )
            {
                double thisGenotypeLogLikelihood = 0.0;
                const auto hapIndicies = genotypes.getHaplotypeIndices( genotypeIndex );

                // prior probability of any haplotype (given genotype).  Assumes uniform distribution across haplotypes
                const auto haplotypePriorProbability = 1.0 / static_cast< double >( hapIndicies.size() );

                for ( std::size_t readIndex = 0; readIndex < nReads; ++readIndex )
                {
                    double sumReadLikelihoodsOverHaplotypes = 0.0;
                    for ( const auto & hapIndex : hapIndicies )
                    {
                        sumReadLikelihoodsOverHaplotypes += probReadsGivenHaplotypes( readIndex, hapIndex );
                    }
                    thisGenotypeLogLikelihood +=
                        std::log( sumReadLikelihoodsOverHaplotypes * haplotypePriorProbability );
                }

                genotypeLikelihoods[genotypeIndex] = thisGenotypeLogLikelihood;
            }

            // Now re-scale all the genotype likelihoods for this individual, to prevent underflows, and
            rescaleLogLikelihoods( genotypeLikelihoods, nReads > 0 );

            // store the likelihoods as non-log values
            std::transform( genotypeLikelihoods.begin(), genotypeLikelihoods.end(), genotypeLikelihoods.begin(), exp );

            if ( false )
            {
                // print likelihoods of genotypes
                ECHIDNA_LOG( SUPER_DEBUG, "Genotype likelihoods:-" );
                for ( std::size_t genotypeIndex = 0; genotypeIndex < genotypes.size(); ++genotypeIndex )
                {
                    if ( genotypeLikelihoods[genotypeIndex] > 1e-6 )
                    {
                        ECHIDNA_LOG( SUPER_DEBUG, std::to_string( genotypeLikelihoods[genotypeIndex] ) + " for " +
                                                      genotypes[genotypeIndex]->toString( haplotypes ) );
                    }
                }
            }

            return genotypeLikelihoods;
        }
    }
}
}
