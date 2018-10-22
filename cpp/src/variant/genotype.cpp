// All content Copyright (C) 2018 Genomics plc
#include <cassert>
#include <algorithm>
#include "variant/genotype.hpp"
#include <unordered_map>
#include "utils/combinationGenerator.hpp"
#include "utils/multinomialCoefficients.hpp"

namespace wecall
{
namespace variant
{

    std::size_t computeNCombinationsThisGenotype( const haplotypeAndCount_t & haplotypes )
    {
        // extract the vector of haplotype counts
        std::vector< unsigned int > counts( haplotypes.size() );
        for ( auto & keyValue : haplotypes )
        {
            counts.push_back( keyValue.second );
        }

        return multinomial_coefficient( counts );
    }

    Genotype::Genotype( const haplotypeAndCount_t & haplotypes )
        : m_haplotypes( haplotypes ), m_nCombinationsThisGenotype( computeNCombinationsThisGenotype( haplotypes ) )
    {
    }

    nonPhasedGenotypeHash_t Genotype::getVariantNonPhasedHash( const HaplotypeVector & haplotypes,
                                                               const std::vector< varPtr_t > & variants ) const
    {
        nonPhasedGenotypeHash_t genotypeHash( variants.size(), 0 );
        for ( std::size_t varIndex = 0; varIndex < variants.size(); ++varIndex )
        {
            genotypeHash[varIndex] = static_cast< char >( nStrandsContainingVariant( haplotypes, variants[varIndex] ) );
        }
        return genotypeHash;
    }

    std::size_t Genotype::nStrandsContainingVariant( const HaplotypeVector & haplotypes,
                                                     const variant::varPtr_t & var ) const
    {
        std::size_t numberOfStrands = 0;
        for ( const auto & haplotype : m_haplotypes )
        {
            if ( haplotypes[haplotype.first].containsVariant( var ) )
            {
                numberOfStrands += haplotype.second;
            }
        }
        return numberOfStrands;
    }

    std::string Genotype::toString() const
    {
        std::stringstream message;
        if ( m_haplotypes.size() == 1 )
        {
            message << "Homozygous [ ";
        }
        else if ( m_haplotypes.size() == 2 )
        {
            message << "Heterozygous [ ";
        }
        else
        {
            message << m_haplotypes.size() << "-ploid [";
        }
        for ( auto hapCount : m_haplotypes )
        {
            message << hapCount.second << " * " << hapCount.first << ", ";
        }
        message << " ] ";
        return message.str();
    }

    std::string Genotype::toString( const HaplotypeVector & haplotypes ) const
    {
        std::stringstream message;
        message << "[ ";
        for ( auto hapCount : m_haplotypes )
        {
            message << hapCount.second << "*" << haplotypes[hapCount.first] << ", ";
        }
        message << " ] ";
        return message.str();
    }

    genoCall_t Genotype::getGenotypeCall( const HaplotypeVector & haplotypes, const variant::varPtr_t & var ) const
    {
        genoCall_t genoCall;
        for ( const auto & haplotype : m_haplotypes )
        {
            std::fill_n( std::back_inserter( genoCall ), haplotype.second,
                         haplotypes[haplotype.first].getCallType( var ) );
        }
        return genoCall;
    }

    std::vector< unsigned int > Genotype::getHaplotypeIndices() const
    {
        std::vector< unsigned int > haplotypeIndizes;
        for ( const auto & haplotype : m_haplotypes )
        {
            std::fill_n( std::back_inserter( haplotypeIndizes ), haplotype.second, haplotype.first );
        }
        return haplotypeIndizes;
    }

    GenotypeVector::GenotypeVector( const std::size_t ploidy,
                                    const HaplotypeVector & haplotypes,
                                    const std::vector< std::size_t > & hapIndiciesToUse,
                                    const std::vector< variant::varPtr_t > & variants )
        : m_ploidy( ploidy ), m_variants( variants ), m_genotypes(), m_genotypeHaplotypeIndexes()
    {
        if ( m_ploidy > 0 )
        {
            for ( std::size_t hapIndexPtr1 = 0; hapIndexPtr1 != hapIndiciesToUse.size(); ++hapIndexPtr1 )
            {
                const auto & hap1 = haplotypes[hapIndiciesToUse[hapIndexPtr1]];
                for ( std::size_t hapIndexPtr2 = hapIndexPtr1 + 1; hapIndexPtr2 != hapIndiciesToUse.size();
                      ++hapIndexPtr2 )
                {
                    const auto & hap2 = haplotypes[hapIndiciesToUse[hapIndexPtr2]];
                    ECHIDNA_ASSERT( hap1 != hap2, "Can not form GenotypeVector with an unmerged HaplotypeVector" );
                }
            }

            for ( const std::vector< unsigned long > & indexPtrs :
                  CombinationGenerator( m_ploidy, hapIndiciesToUse.size() ) )
            {
                assert( ploidy == indexPtrs.size() );
                std::vector< unsigned long > indices( ploidy );

                std::unordered_map< unsigned long, unsigned long > indiciesToCounts;

                for ( std::size_t i = 0; i < indexPtrs.size(); ++i )
                {
                    const auto newIndex = hapIndiciesToUse[indexPtrs[i]];
                    indices[i] = newIndex;

                    if ( indiciesToCounts.count( newIndex ) )
                    {
                        ++indiciesToCounts[newIndex];
                    }
                    else
                    {
                        indiciesToCounts[newIndex] = 1;
                    }
                }

                haplotypeAndCount_t hapPointers;
                for ( const auto & idx : indiciesToCounts )
                {
                    hapPointers.push_back( std::make_pair( idx.first, idx.second ) );
                }

                m_genotypes.emplace_back( std::make_shared< variant::Genotype >( hapPointers ) );
                m_genotypeHaplotypeIndexes.push_back( indices );
            }
            this->buildGenotypeWithSameNonPhaseRepresentation( haplotypes );
        }
    }

    std::set< std::size_t > GenotypeVector::genotypesWithSameNonPhasedRepresentation(
        const std::size_t givenGenotypeIndex ) const
    {
        ECHIDNA_ASSERT( givenGenotypeIndex < m_genotypes.size(), "Genotype index out-of-range" );
        ECHIDNA_ASSERT( m_genotypeWithSameNonPhasedRepresentation.find( givenGenotypeIndex ) !=
                            m_genotypeWithSameNonPhasedRepresentation.end(),
                        "Genotype index not found in non-phased representation map" );

        return *( m_genotypeWithSameNonPhasedRepresentation.at( givenGenotypeIndex ) );
    }

    void GenotypeVector::buildGenotypeWithSameNonPhaseRepresentation( const HaplotypeVector & haplotypes )
    {
        std::unordered_map< nonPhasedGenotypeHash_t, genotypeIndexSetPtr_t > nonPhasedGenotypeSets;

        for ( std::size_t i = 0; i < m_genotypes.size(); ++i )
        {
            auto const hash = m_genotypes[i]->getVariantNonPhasedHash( haplotypes, m_variants );

            if ( nonPhasedGenotypeSets.find( hash ) == nonPhasedGenotypeSets.end() )
            {
                nonPhasedGenotypeSets[hash] = std::make_shared< genotypeIndexSet_t >();
            }
            nonPhasedGenotypeSets[hash]->insert( i );
        }

        for ( const auto & hashAndSet : nonPhasedGenotypeSets )
        {
            for ( const auto genotypeIndex : *( hashAndSet.second ) )
            {
                m_genotypeWithSameNonPhasedRepresentation[genotypeIndex] = hashAndSet.second;
            }
        }
    }
}
}
