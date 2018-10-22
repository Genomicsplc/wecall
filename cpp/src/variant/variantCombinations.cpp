// All content Copyright (C) 2018 Genomics plc
#include "haplotype.hpp"
#include "variant/clustering.hpp"
#include "variant/variantCombinations.hpp"
#include "io/read.hpp"

namespace wecall
{
namespace variant
{
    VariantCombinations::VariantCombinations( const std::size_t minReadsToSupportClaim,
                                              const int64_t maxClusterDistance )
        : m_variantCombinations(),
          m_minReadsToSupportClaim( minReadsToSupportClaim ),
          m_maxClusterDistance( maxClusterDistance )
    {
    }

    void VariantCombinations::computeVariantCombinations( const std::vector< varPtr_t > & variants,
                                                          const std::size_t maxCombinations,
                                                          const utils::referenceSequencePtr_t & reference )
    {
        m_allCombinationsComputed = true;

        // if there are no variants, only add the reference combo
        if ( variants.empty() )
        {
            m_variantCombinations.push_back( {} );
            return;
        }

        // init variant combinations with first variant
        auto it = variants.begin();
        auto first_var = *it;
        m_variantCombinations.push_back( {first_var} );
        it++;

        // add valid combinations for subsequent variants in cluster
        for ( ; it != variants.end(); it++ )
        {
            // include the reference sequence when checking for max number of combinations
            if ( m_variantCombinations.size() + 1 >= maxCombinations )
            {
                m_allCombinationsComputed = false;
                m_variantCombinations.clear();
                ECHIDNA_LOG( DEBUG, "Too many haplotype combinations." );
                return;
            }

            auto variant = *it;
            std::vector< std::vector< varPtr_t > > newCombos;
            std::vector< varPtr_t > alwaysTogetherVariants;
            std::vector< varPtr_t > neverTogetherVariants;
            std::vector< varPtr_t > secondImpliesFirstVariants;
            for ( auto & prevCombo : m_variantCombinations )
            {
                if ( isValidCombinationVec( {prevCombo.back(), variant}, reference ) )
                {
                    VariantPairCombinationState is_related = getState( prevCombo.back(), variant );
                    if ( is_related == VariantPairCombinationState::ALWAYS_TOGETHER )
                    {
                        // append variant to current combination
                        alwaysTogetherVariants.emplace_back( prevCombo.back() );
                        prevCombo.push_back( variant );
                    }
                    else if ( is_related == VariantPairCombinationState::NEVER_TOGETHER )
                    {
                        neverTogetherVariants.emplace_back( prevCombo.back() );
                    }
                    else if ( is_related == VariantPairCombinationState::FIRST_IMPLIES_SECOND )
                    {
                        // append variant to current combination
                        prevCombo.push_back( variant );
                    }
                    else if ( is_related == VariantPairCombinationState::SECOND_IMPLIES_FIRST )
                    {
                        // add the combination of prevCombo + variant
                        secondImpliesFirstVariants.emplace_back( prevCombo.back() );
                        auto newCombo = prevCombo;
                        newCombo.push_back( variant );
                        newCombos.push_back( newCombo );
                    }
                    else
                    {
                        // add the combination of prevCombo + variant
                        auto newCombo = prevCombo;
                        newCombo.push_back( variant );
                        newCombos.push_back( newCombo );
                    }
                }
            }
            // add the current variant once (will be filtered if it violates the never together case)
            if ( alwaysTogetherVariants.empty() and secondImpliesFirstVariants.empty() )
            {
                newCombos.push_back( {variant} );
            }

            auto validNewCombos =
                filterVariantCombinations( newCombos, alwaysTogetherVariants, neverTogetherVariants, variant );
            m_variantCombinations.insert( m_variantCombinations.end(), validNewCombos.begin(), validNewCombos.end() );
        }
        m_variantCombinations.push_back( {} );  // add reference combo.

        ECHIDNA_LOG( DEBUG, "Generated " << m_variantCombinations.size() << " haplotypes." );
        if ( false )
        {
            std::stringstream debugMessage;
            debugMessage << "\t";
            for ( auto & combo : m_variantCombinations )
            {
                for ( auto & var : combo )
                {
                    debugMessage << var->toString();
                    debugMessage << ", ";
                }
                debugMessage << std::endl
                             << "\t";
            }
            ECHIDNA_LOG( SUPER_DEBUG, debugMessage.str() );
        }
    }

    std::vector< std::vector< varPtr_t > > VariantCombinations::filterVariantCombinations(
        const std::vector< std::vector< varPtr_t > > & variantCombinations,
        const std::vector< varPtr_t > & alwaysTogetherVariants,
        const std::vector< varPtr_t > & neverTogetherVariants,
        varPtr_t currentVariant )
    {
        // check that there are no invalid combinations due to always together in the new Combos
        if ( alwaysTogetherVariants.empty() and neverTogetherVariants.empty() )
        {
            return variantCombinations;
        }

        std::vector< std::vector< varPtr_t > > validVariantCombinations;
        for ( auto & combo : variantCombinations )
        {
            bool invalid = false;
            bool currentVariantExists = std::find( combo.begin(), combo.end(), currentVariant ) != combo.end();

            // check that it is valid against always together variants
            for ( auto & variant : alwaysTogetherVariants )
            {
                bool relatedVarExists = std::find( combo.begin(), combo.end(), variant ) != combo.end();
                if ( ( currentVariantExists and ( not relatedVarExists ) ) or
                     ( ( not currentVariantExists ) and relatedVarExists ) )
                {
                    invalid = true;
                }
            }

            for ( auto & variant : neverTogetherVariants )
            {
                bool relatedVarExists = std::find( combo.begin(), combo.end(), variant ) != combo.end();
                if ( currentVariantExists and relatedVarExists )
                {
                    invalid = true;
                }
            }

            if ( not invalid )
            {
                validVariantCombinations.push_back( combo );
            }
        }
        return validVariantCombinations;
    }

    VariantPairCombinationState VariantCombinations::getState( varPtr_t first, varPtr_t second ) const
    {
        auto firstInterval = first->getStartEndRegions().getSpan();
        auto secondInterval = second->getStartEndRegions().getSpan();
        firstInterval.combine( secondInterval );
        if ( static_cast< std::size_t >( firstInterval.size() ) >= m_maxClusterDistance )
        {
            return VariantPairCombinationState::UNCERTAIN;
        }

        std::vector< io::readPtr_t > firstReads;
        std::vector< io::readPtr_t > secondReads;
        getOverlappingReads( first, second, firstReads, secondReads );

        // find intersection between reads of the first and second variant
        std::sort( firstReads.begin(), firstReads.end() );
        std::sort( secondReads.begin(), secondReads.end() );

        std::vector< io::readPtr_t > shared_reads( std::max( firstReads.size(), secondReads.size() ) );
        std::vector< io::readPtr_t >::iterator it;

        it = std::set_intersection( firstReads.begin(), firstReads.end(), secondReads.begin(), secondReads.end(),
                                    shared_reads.begin() );
        shared_reads.resize( it - shared_reads.begin() );

        // compute number of reads supporting different scenarios
        std::size_t numAlwaysTogether = shared_reads.size();
        std::size_t numOnlyVar1 = firstReads.size();
        std::size_t numOnlyVar2 = secondReads.size();

        if ( numAlwaysTogether == 0 and numOnlyVar1 >= m_minReadsToSupportClaim and
             numOnlyVar2 >= m_minReadsToSupportClaim )
        {
            return VariantPairCombinationState::NEVER_TOGETHER;
        }
        if ( numAlwaysTogether == numOnlyVar1 and numAlwaysTogether == numOnlyVar2 and
             numAlwaysTogether >= m_minReadsToSupportClaim )
        {
            return VariantPairCombinationState::ALWAYS_TOGETHER;
        }
        if ( numAlwaysTogether == numOnlyVar2 and numOnlyVar1 > numOnlyVar2 and
             numOnlyVar2 >= m_minReadsToSupportClaim )
        {
            return VariantPairCombinationState::SECOND_IMPLIES_FIRST;
        }
        if ( numAlwaysTogether == numOnlyVar1 and numOnlyVar2 > numOnlyVar1 and
             numOnlyVar1 >= m_minReadsToSupportClaim )
        {
            return VariantPairCombinationState::FIRST_IMPLIES_SECOND;
        }
        return VariantPairCombinationState::UNCERTAIN;
    }

    void VariantCombinations::getOverlappingReads( const varPtr_t & first,
                                                   const varPtr_t & second,
                                                   std::vector< io::readPtr_t > & firstReads,
                                                   std::vector< io::readPtr_t > & secondReads ) const
    {
        firstReads = first->getReads();
        secondReads = second->getReads();
        const auto firstInterval = first->getStartEndRegions().getSpan().interval();
        const auto secondInterval = second->getStartEndRegions().getSpan().interval();

        auto filterFunc = [firstInterval, secondInterval]( io::readPtr_t readPtr )
        {
            const auto readInterval = readPtr->getMaximalReadInterval();
            return not( readInterval.contains( firstInterval ) and readInterval.contains( secondInterval ) );
        };

        firstReads.erase( remove_if( firstReads.begin(), firstReads.end(), filterFunc ), firstReads.end() );
        secondReads.erase( remove_if( secondReads.begin(), secondReads.end(), filterFunc ), secondReads.end() );
    }
}
}