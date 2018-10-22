// All content Copyright (C) 2018 Genomics plc
#include "variant/haplotype.hpp"
#include "caller/region.hpp"
#include "utils/interval.hpp"
#include "utils/logging.hpp"
#include "utils/exceptions.hpp"
#include <algorithm>
#include "utils/partition.hpp"
#include "variant/type/variant.hpp"
#include "variantNormalizer.hpp"

namespace wecall
{
namespace variant
{

    bool isValidCombinationVec( const std::vector< varPtr_t > & variants,
                                utils::referenceSequencePtr_t referenceSequence )
    {
        return isValidVariantCombination( variants.cbegin(), variants.cend(), referenceSequence );
    }

    double HaplotypeScore::score() const
    {
        return this->readKmersExplainedByHaploltype * ( 1.0 - this->kmersUniqueToHaplotype ) +
               minReadDataCountToHaplotypeCount;
    }

    Haplotype::Haplotype( utils::referenceSequencePtr_t referenceSequence,
                          const caller::SetRegions & hapRegions,
                          const variantSet_t & variants,
                          const std::size_t referencePadding )
        : m_referenceSequence( referenceSequence ),
          m_regions( hapRegions ),
          m_paddedSequences(),
          m_variants( variants ),
          m_variantsPerRegion(),
          m_referencePadding( referencePadding )
    {
        ECHIDNA_ASSERT( isValidVariantCombination( variants.cbegin(), variants.cend(), referenceSequence ),
                        "Attempted to form invalid haplotype: " + this->toString() );

        for ( const auto & var : m_variants )
        {
            m_regions.insert( var->region() );
        }

        m_regions.fill( m_referencePadding );

        auto currentVariant = variants.begin();
        for ( const auto & miniHaplotypeRegion : m_regions )
        {
            variantSet_t current;
            while ( currentVariant != variants.end() and miniHaplotypeRegion.contains( ( *currentVariant )->region() ) )
            {
                current.insert( *currentVariant );
                ++currentVariant;
            }
            const caller::Region left( m_referenceSequence->contig(), miniHaplotypeRegion.start() - m_referencePadding,
                                       miniHaplotypeRegion.start() );
            const caller::Region right( m_referenceSequence->contig(), miniHaplotypeRegion.end(),
                                        miniHaplotypeRegion.end() + m_referencePadding );

            m_variantsPerRegion.emplace_back( current );
            m_paddedSequences.emplace_back( m_referenceSequence->subseq( left ).sequence() +
                                            buildHaplotypeSequence( referenceSequence, miniHaplotypeRegion, current ) +
                                            m_referenceSequence->subseq( right ).sequence() );
        }

        if ( m_paddedSequences.front().size() > 1000 )
        {
            ECHIDNA_LOG( DEBUG, "Large haplotype constructed of size " << m_paddedSequences.front().size() );
            ECHIDNA_LOG( DEBUG, "m_startPos = " << m_referenceSequence->start()
                                                << " and m_endPos = " << m_referenceSequence->end() );
        }
    }

    varPtr_t mnpFromSNPs( utils::referenceSequencePtr_t reference, const std::vector< varPtr_t > & snps )
    {
        const caller::Region region( reference->contig(), snps.front()->start(), snps.back()->end() );
        const auto alt =
            Haplotype::buildHaplotypeSequence( reference, region, variantSet_t( snps.cbegin(), snps.cend() ) );
        return std::make_shared< Variant >( reference, region, alt, true );
    }

    variantSet_t Haplotype::withNormalizedVariants()
    {
        VariantNormalizer variantNormalizer( m_referenceSequence );

        variantSet_t allVariants;
        std::vector< variantSet_t > variantsPerRegion;

        assert( m_regions.size() == m_paddedSequences.size() );
        assert( m_regions.size() == m_variantsPerRegion.size() );

        std::size_t index = 0;
        for ( const auto & region : m_regions )
        {
            const auto & paddedSequence = m_paddedSequences[index];
            const auto & variantsThisRegion = m_variantsPerRegion[index];

            const auto paddedRegion = region.getPadded( m_referencePadding );
            const auto normalVariants =
                variantNormalizer.getNormalized( paddedRegion, paddedSequence, variantsThisRegion );

            variantsPerRegion.emplace_back( normalVariants.cbegin(), normalVariants.cend() );
            allVariants.insert( normalVariants.cbegin(), normalVariants.cend() );
            ++index;
        }
        m_variantsPerRegion = variantsPerRegion;
        m_variants = allVariants;
        return m_variants;
    }

    variantSet_t Haplotype::withMNPs()
    {
        const auto indelComp = []( varPtr_t var1, varPtr_t var2 )
        {
            return var1->isIndel() == var2->isIndel();
        };

        const auto setVariants = this->getVariants();
        const std::vector< varPtr_t > allVariants( setVariants.begin(), setVariants.end() );

        const auto partitionedVariants = utils::functional::partition( allVariants, indelComp );

        m_variants.clear();
        variantSet_t newMNPs;
        for ( const auto & partition : partitionedVariants )
        {
            if ( partition.front()->isIndel() )
            {
                m_variants.insert( partition.begin(), partition.end() );
            }
            else if ( partition.size() == 1 )
            {
                m_variants.insert( partition.front() );
            }
            else
            {
                const auto newMNP = mnpFromSNPs( this->m_referenceSequence, partition );
                m_variants.insert( newMNP );
                newMNPs.insert( newMNP );
            }
        }
        return newMNPs;
    }

    utils::BasePairSequence Haplotype::buildHaplotypeSequence( const utils::referenceSequencePtr_t & referenceSequence,
                                                               const caller::Region haplotypeRegion,
                                                               const variantSet_t & variants )
    {
        int64_t currentPos = haplotypeRegion.start();

        std::string hapSequence;
        for ( auto varPtr : variants )
        {
            const caller::Region regionBeforeVariant( haplotypeRegion.contig(), currentPos, varPtr->start() );
            hapSequence += referenceSequence->subseq( regionBeforeVariant ).sequence().str();
            hapSequence += varPtr->sequence().str();
            currentPos = varPtr->end();
        }
        const caller::Region regionAfterVariants( haplotypeRegion.contig(), currentPos, haplotypeRegion.end() );
        hapSequence += referenceSequence->subseq( regionAfterVariants ).sequence().str();

        return hapSequence;
    }

    //-----------------------------------------------------------------------------------------

    Haplotype::~Haplotype() {}

    //-----------------------------------------------------------------------------------------

    std::string Haplotype::toString() const
    {
        std::stringstream message;
        message << "Haplotype(";

        for ( const auto & region : m_regions )
        {
            message << region.toString() << ", ";
        }
        message << "{";

        for ( auto varPtr : m_variants )
        {
            message << varPtr->toString() << ", ";
        }

        message << "})";
        return message.str();
    }

    //-----------------------------------------------------------------------------------------

    double Haplotype::productOfVariantPriors() const
    {
        auto product_of_priors = 1.0;

        for ( auto var_ptr : m_variants )
        {
            product_of_priors *= var_ptr->prior();
        }
        return product_of_priors;
    }

    bool Haplotype::isReference( const caller::Region & region ) const
    {
        const auto overlaps = [&region]( varPtr_t varPtr ) -> bool
        {
            const auto varOverlaps = varPtr->overlaps( region );
            if ( not varOverlaps )
            {
                return false;
            }
            else if ( varPtr->isIndel() or varPtr->sequenceLengthInRef() == 1 )
            {
                return varOverlaps;
            }
            else
            {
                const auto start = std::max( region.start(), varPtr->start() );
                const auto size = std::min( region.end(), varPtr->end() ) - start;

                const auto varRef = varPtr->refSequence().sequence().substr( start - varPtr->start(), size );
                const auto varMut = varPtr->sequence().substr( start - varPtr->start(), size );
                return varRef != varMut;
            }
        };

        return not std::any_of( m_variants.begin(), m_variants.end(), overlaps );
    }

    bool Haplotype::operator<( const Haplotype & other ) const
    {
        if ( m_regions.size() != other.m_regions.size() )
        {
            return m_regions.size() < other.m_regions.size();
        }
        else if ( m_regions != other.m_regions )
        {
            return m_regions < other.m_regions;
        }
        else if ( paddedSequences().size() != other.paddedSequences().size() )
        {
            return ( paddedSequences().size() < other.paddedSequences().size() );
        }

        for ( std::size_t i = 0; i < paddedSequences().size(); ++i )
        {
            if ( paddedSequences()[i] != other.paddedSequences()[i] )
            {
                return paddedSequences()[i] < other.paddedSequences()[i];
            }
        }
        return false;
    }

    bool Haplotype::operator==( const Haplotype & rhs ) const
    {
        return m_regions == rhs.m_regions and m_paddedSequences.size() == rhs.m_paddedSequences.size() and
               std::equal( m_paddedSequences.begin(), m_paddedSequences.end(), rhs.m_paddedSequences.begin() );
    }

    caller::Call::Type Haplotype::getCallType( const varPtr_t & varPtr ) const
    {
        if ( this->containsVariant( varPtr ) )
        {
            return caller::Call::VAR;
        }
        else if ( this->isReference( varPtr->region() ) )
        {
            return caller::Call::REF;
        }
        else
        {
            return caller::Call::UNKNOWN;
        }
    }

    //-----------------------------------------------------------------------------------------

    bool Haplotype::isMoreLikelyThan( const Haplotype & rhs ) const
    {
        const auto fromBreakpoint = []( const varPtr_t & varPtr ) -> bool
        {
            return varPtr->isFromBreakpoint();
        };
        const auto lhsVariants = this->getVariants();
        const auto rhsVariants = rhs.getVariants();
        const bool lhsFromBp = std::any_of( lhsVariants.cbegin(), lhsVariants.cend(), fromBreakpoint );
        const bool rhsFromBp = std::any_of( rhsVariants.cbegin(), rhsVariants.cend(), fromBreakpoint );
        if ( lhsFromBp != rhsFromBp )
        {
            return not lhsFromBp;
        }

        const auto priorDiff = this->productOfVariantPriors() - rhs.productOfVariantPriors();

        if ( std::abs( priorDiff ) > std::numeric_limits< double >::epsilon() )
        {
            return priorDiff > 0;
        }

        else if ( this->numberOfVariants() != rhs.numberOfVariants() )
        {
            return this->numberOfVariants() < rhs.numberOfVariants();
        }

        else if ( this->paddedSequences().size() != rhs.paddedSequences().size() )
        {
            return this->paddedSequences().size() < rhs.paddedSequences().size();
        }
        else
        {
            for ( std::size_t i = 0; i < this->paddedSequences().size(); ++i )
            {
                if ( this->paddedSequences()[i] < rhs.paddedSequences()[i] )
                {
                    return true;
                }
            }
            return false;
        }
    }

    HaplotypeVector::HaplotypeVector( caller::SetRegions regions, utils::referenceSequencePtr_t referenceSequence )
        : m_regions( regions ), m_paddedReferenceSequence( referenceSequence )
    {
        ECHIDNA_ASSERT( m_paddedReferenceSequence->region().contains( m_regions.getSpan() ),
                        "Require compatible reference" );
        const auto span = m_regions.getSpan();
        m_referencePadding = int64_to_sizet(
            std::min( span.start() - referenceSequence->start(), referenceSequence->end() - span.end() ) );
        m_regions.fill( m_referencePadding );
    }

    //-----------------------------------------------------------------------------------------

    std::set< std::size_t > HaplotypeVector::getHaplotypeIndicesForVariant( const varPtr_t & variant ) const
    {
        std::set< std::size_t > hapIndicesForThisVariant;

        for ( std::size_t hapIndex = 0; hapIndex < m_haplotypes.size(); ++hapIndex )
        {
            if ( m_haplotypes[hapIndex].containsVariant( variant ) )
            {
                hapIndicesForThisVariant.insert( hapIndex );
            }
        }

        return hapIndicesForThisVariant;
    }

    std::set< std::size_t > HaplotypeVector::getIndicesForReference( const caller::Region & region ) const
    {
        std::set< std::size_t > hapIndicesForReferenceAtRegion;

        for ( std::size_t hapIndex = 0; hapIndex < m_haplotypes.size(); ++hapIndex )
        {
            if ( m_haplotypes[hapIndex].isReference( region ) )
            {
                hapIndicesForReferenceAtRegion.insert( hapIndex );
            }
        }
        return hapIndicesForReferenceAtRegion;
    }

    //-----------------------------------------------------------------------------------------

    void HaplotypeVector::sort() { std::sort( m_haplotypes.begin(), m_haplotypes.end() ); }

    std::string HaplotypeVector::toString() const
    {
        std::stringstream s;
        for ( const auto & haplotype : m_haplotypes )
        {
            s << haplotype.toString() << std::endl;
        }
        return s.str();
    }

    void HaplotypeVector::merge()
    {
        //        m_haplotypes.erase( std::unique( m_haplotypes.begin(), m_haplotypes.end() ), m_haplotypes.end() );

        std::vector< Haplotype > previousHaplotypes = m_haplotypes;

        m_haplotypes.clear();

        for ( auto & haplotype : previousHaplotypes )
        {
            if ( m_haplotypes.empty() )
            {
                m_haplotypes.push_back( haplotype );
            }

            else
            {
                const auto prevHaplotype = m_haplotypes.back();

                if ( prevHaplotype != haplotype )
                {
                    m_haplotypes.push_back( haplotype );
                }

                else
                {
                    const bool keepHap = haplotype.isMoreLikelyThan( prevHaplotype );
                    m_haplotypes.back() = keepHap ? haplotype : prevHaplotype;
                }
            }
        }
    }

    variantSet_t HaplotypeVector::withMNPs()
    {
        variantSet_t newMNPs;
        for ( auto & haplotype : m_haplotypes )
        {
            const auto hapNewMNPS = haplotype.withMNPs();
            newMNPs.insert( hapNewMNPS.cbegin(), hapNewMNPS.cend() );
        }
        return newMNPs;
    }

    variantSet_t HaplotypeVector::withNormalizedVariants()
    {
        variantSet_t normalVariants;
        for ( auto & haplotype : m_haplotypes )
        {
            const auto hapNewVariants = haplotype.withNormalizedVariants();
            normalVariants.insert( hapNewVariants.cbegin(), hapNewVariants.cend() );
        }
        return normalVariants;
    }

    void HaplotypeVector::push_back( const variantSet_t & varCombo )
    {
        m_haplotypes.emplace_back( m_paddedReferenceSequence, m_regions, varCombo, m_referencePadding );
    }

    void HaplotypeVector::push_back( const variantSet_t & varCombo, const size_t haplotypeId )
    {
        Haplotype haplotype = Haplotype( m_paddedReferenceSequence, m_regions, varCombo, m_referencePadding );
        haplotype.setId( haplotypeId );
        m_haplotypes.push_back( haplotype );
    }

    void HaplotypeVector::keepIndicies( std::set< std::size_t > indicies )
    {
        std::vector< Haplotype > retained;

        for ( const auto & keep : indicies )
        {
            retained.push_back( this->m_haplotypes[keep] );
        }
        m_haplotypes = retained;
    }
}
}
