// All content Copyright (C) 2018 Genomics plc
#ifndef GENOTYPE_HPP
#define GENOTYPE_HPP

#include <cstdlib>
#include <memory>

#include <unordered_map>

#include "variant/haplotype.hpp"
#include "caller/callSet.hpp"

namespace echidna
{
namespace variant
{
    using haplotypeAndCount_t = std::vector< std::pair< std::size_t, unsigned int > >;
    using nonPhasedGenotypeHash_t = std::string;

    class Genotype
    {
    public:
        Genotype( const haplotypeAndCount_t & haplotypes );

        std::size_t nStrandsContainingVariant( const HaplotypeVector & haplotypes,
                                               const variant::varPtr_t & var ) const;
        std::size_t nCombinationsThisGenotype() const { return m_nCombinationsThisGenotype; }
        std::string toString() const;
        std::string toString( const HaplotypeVector & ) const;
        genoCall_t getGenotypeCall( const HaplotypeVector & haplotypes, const variant::varPtr_t & var ) const;
        std::vector< unsigned int > getHaplotypeIndices() const;

        nonPhasedGenotypeHash_t getVariantNonPhasedHash( const HaplotypeVector & haplotypes,
                                                         const std::vector< varPtr_t > & variants ) const;

    private:
        const haplotypeAndCount_t m_haplotypes;
        const std::size_t m_nCombinationsThisGenotype;
    };

    typedef std::shared_ptr< Genotype > genotypePtr_t;

    class GenotypeVector
    {
    public:
        using genotypeIndexSet_t = std::set< std::size_t >;
        using genotypeIndexSetPtr_t = std::shared_ptr< genotypeIndexSet_t >;
        using genotypeWithSameNonPhasedRepresentation_t = std::unordered_map< std::size_t, genotypeIndexSetPtr_t >;

    public:
        GenotypeVector( const std::size_t ploidy,
                        const HaplotypeVector & haplotypes,
                        const std::vector< std::size_t > & hapIndiciesToUse,
                        const std::vector< variant::varPtr_t > & variants );

        std::vector< std::size_t > getHaplotypeIndices( const std::size_t genotypeIndex ) const
        {
            return m_genotypeHaplotypeIndexes.at( genotypeIndex );
        }

        std::set< std::size_t > genotypesWithSameNonPhasedRepresentation( const std::size_t givenGenotypeIndex ) const;

        std::size_t size() const { return m_genotypes.size(); }
        std::size_t ploidy() const { return m_ploidy; }

        genotypePtr_t operator[]( std::size_t genotypeIndex ) const { return m_genotypes.at( genotypeIndex ); }

    private:
        void buildGenotypeWithSameNonPhaseRepresentation( const HaplotypeVector & haplotypes );

    private:
        const std::size_t m_ploidy;
        const std::vector< varPtr_t > m_variants;

        std::vector< genotypePtr_t > m_genotypes;

        std::vector< std::vector< std::size_t > > m_genotypeHaplotypeIndexes;

        genotypeWithSameNonPhasedRepresentation_t m_genotypeWithSameNonPhasedRepresentation;
    };
}
}

#endif
