// All content Copyright (C) 2018 Genomics plc
#ifndef HAPLOTYPE_HPP
#define HAPLOTYPE_HPP

#include <caller/callSet.hpp>
#include "common.hpp"

#include "variant/clustering.hpp"
#include "variant/type/variant.hpp"
#include "mapping/hashMapper.hpp"
#include "io/fastaFile.hpp"
#include "utils/logging.hpp"
#include "utils/interval.hpp"

namespace echidna
{
namespace variant
{
    struct HaplotypeScore
    {
        double readKmersExplainedByHaploltype;
        double kmersUniqueToHaplotype;
        double minReadDataCountToHaplotypeCount;

        double score() const;
    };

    using hapScorePtr_t = std::shared_ptr< HaplotypeScore >;

    class Haplotype
    {
    public:
        Haplotype( utils::referenceSequencePtr_t referenceSequence,
                   const caller::SetRegions & hapRegions,
                   const variantSet_t & variants,
                   const std::size_t referencePadding );

        ~Haplotype();

        const caller::SetRegions & regions() const { return m_regions; }
        const std::vector< utils::BasePairSequence > & paddedSequences() const { return m_paddedSequences; }

        bool operator<( const Haplotype & other ) const;
        bool operator==( const Haplotype & rhs ) const;
        bool operator!=( const Haplotype & rhs ) const { return not( *this == rhs ); }

        bool containsVariant( const varPtr_t & varPtr ) const { return ( m_variants.count( varPtr ) > 0 ); }
        caller::Call::Type getCallType( const varPtr_t & varPtr ) const;

        bool isReference( const caller::Region & region ) const;

        bool isMoreLikelyThan( const Haplotype & rhs ) const;
        double productOfVariantPriors() const;
        int numberOfVariants() const { return static_cast< int >( m_variants.size() ); }
        std::string toString() const;

        friend std::ostream & operator<<( std::ostream & out, const Haplotype & haplotype )
        {
            return out << haplotype.toString();
        }

        variantSet_t withMNPs();
        variantSet_t withNormalizedVariants();

        variantSet_t getVariants() const { return m_variants; }

        size_t getId() const { return m_id; }
        void setId( const size_t id ) { m_id = id; }

        static utils::BasePairSequence buildHaplotypeSequence( const utils::referenceSequencePtr_t & referenceSequence,
                                                               const caller::Region haplotypeRegion,
                                                               const variantSet_t & variants );

    private:
        utils::referenceSequencePtr_t m_referenceSequence;
        caller::SetRegions m_regions;
        std::vector< utils::BasePairSequence > m_paddedSequences;

        variantSet_t m_variants;
        std::vector< variantSet_t > m_variantsPerRegion;
        std::size_t m_referencePadding;
        std::size_t m_id = -1;
    };

    class HaplotypeVector
    {

    public:
        using iterator = std::vector< Haplotype >::iterator;
        using const_iterator = std::vector< Haplotype >::const_iterator;

    public:
        HaplotypeVector( caller::SetRegions regions, utils::referenceSequencePtr_t referenceSequence );

        const Haplotype & operator[]( std::size_t hapIndex ) const { return m_haplotypes.at( hapIndex ); }

        std::size_t size() const { return m_haplotypes.size(); }

        std::set< std::size_t > getHaplotypeIndicesForVariant( const varPtr_t & variant ) const;
        std::set< std::size_t > getIndicesForReference( const caller::Region & region ) const;

        void push_back( const variantSet_t & varCombo );
        void push_back( const variantSet_t & varCombo, const size_t haplotypeId );
        void keepIndicies( std::set< std::size_t > indicies );
        bool empty() const { return m_haplotypes.empty(); }

        const Haplotype & front() const { return m_haplotypes.front(); }
        const Haplotype & back() const { return m_haplotypes.back(); }

        const_iterator begin() const { return m_haplotypes.cbegin(); }
        const_iterator end() const { return m_haplotypes.cend(); }

        void sort();
        void merge();
        variantSet_t withMNPs();
        variantSet_t withNormalizedVariants();

        std::string toString() const;

        caller::Region region() const { return m_regions.getSpan(); }

        utils::referenceSequencePtr_t paddedReferenceSequence() const { return m_paddedReferenceSequence; }

        utils::referenceSequencePtr_t regionReferenceSequence() const
        {
            return std::make_shared< utils::ReferenceSequence >( m_paddedReferenceSequence->subseq( region() ) );
        }
        utils::referenceSequencePtr_t leftReferencePadding() const
        {
            const caller::Region left( region().contig(), m_paddedReferenceSequence->start(), region().start() );
            return std::make_shared< utils::ReferenceSequence >( m_paddedReferenceSequence->subseq( left ) );
        }
        utils::referenceSequencePtr_t rightReferencePadding() const
        {
            const caller::Region right( region().contig(), region().end(), m_paddedReferenceSequence->end() );
            return std::make_shared< utils::ReferenceSequence >( m_paddedReferenceSequence->subseq( right ) );
        }

    private:
        caller::SetRegions m_regions;
        utils::referenceSequencePtr_t m_paddedReferenceSequence;

        std::vector< Haplotype > m_haplotypes;
        std::size_t m_referencePadding;
    };

    template < typename variantIterator >
    bool isValidVariantCombination( variantIterator varItBegin,
                                    variantIterator varItEnd,
                                    utils::referenceSequencePtr_t referenceSequence )
    {
        if ( std::distance( varItBegin, varItEnd ) <= 1 )
        {
            return true;
        }
        else
        {
            for ( auto it1 = varItBegin, it2 = std::next( it1 ); it2 != varItEnd; ++it1, ++it2 )
            {
                const varPtr_t var1 = ( *it1 );
                const varPtr_t var2 = ( *it2 );

                if ( var1->overlaps( var2 ) )
                {
                    return false;
                }

                if ( not var2->isFullyLeftAligned() )
                {
                    auto refRegion = var1->region();
                    refRegion.combine( var2->region() );
                    const auto hapForVariant1 =
                        Haplotype::buildHaplotypeSequence( referenceSequence, refRegion, {var1} );
                    const auto hapForVariant2 =
                        Haplotype::buildHaplotypeSequence( referenceSequence, refRegion, {var2} );

                    if ( hapForVariant1 == hapForVariant2 )
                    {
                        return false;
                    }
                }
            }

            return true;
        }
    }

    bool isValidCombinationVec( const std::vector< varPtr_t > & variants,
                                utils::referenceSequencePtr_t referenceSequence );
}
}

#endif
