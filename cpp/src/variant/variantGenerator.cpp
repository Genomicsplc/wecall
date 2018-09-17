// All content Copyright (C) 2018 Genomics plc
#include "variant/variantGenerator.hpp"
#include "variant/variantContainer.hpp"

#include <algorithm>
#include <iostream>
#include "utils/partition.hpp"

#include "common.hpp"
#include <deque>
#include <stack>
#include "io/fastaFile.hpp"
#include "variant/type/variant.hpp"

#include "utils/exceptions.hpp"

namespace echidna
{
namespace variant
{
    std::vector< phred_t > getVariantsReadBaseQualities( int64_t startPos,
                                                         const utils::QualitySequence & qualityString,
                                                         const std::vector< varPtr_t > & variants,
                                                         const std::vector< variant::breakpointPtr_t > & breakpoints )
    {
        std::size_t indentIntoRead = 0;
        auto currentRefPosition = startPos;
        std::vector< phred_t > results;
        // Perhaps refactor this!!
        if ( not breakpoints.empty() and breakpoints.front()->pos() == startPos )
        {
            indentIntoRead += breakpoints.front()->sequence().size();
        }

        for ( const auto & var : variants )
        {
            indentIntoRead += var->start() - currentRefPosition;
            ECHIDNA_ASSERT( var->start() >= currentRefPosition, "Non sequential variants added." );
            auto endIndentIntoRead = indentIntoRead + var->sequenceLength();

            if ( indentIntoRead == endIndentIntoRead )
            {
                results.push_back( 1000.0 );
            }
            else
            {
                auto max = *std::max_element( qualityString.begin() + indentIntoRead,
                                              qualityString.begin() + endIndentIntoRead );
                results.push_back( static_cast< phred_t >( max ) );
            }

            indentIntoRead = endIndentIntoRead;
            currentRefPosition = var->end();
        }

        return results;
    }

    std::vector< varPtr_t > normaliseVariantsOnStrand( const std::vector< varPtr_t > & variants,
                                                       const utils::ReferenceSequence & referenceSequence )
    {
        // Mantra: All that is indel must be left-aligned. All that is indel must not have variant on adjacent-left.
        // Left-align takes 1 * indel --> 1 * indel.
        // Join + Split only can produce one indel at maximum. So #outputIndels <= #inputIndels.
        // Left-align and Join don't change total number of ref + alt base-pairs.
        // Split can reduce the total number of ref and alt base-pairs (as split of MNP yields fewer some refs which are
        // skipped)
        // So Sum(#Ref-bps) and Sum(#Alt-bps) is reduced. This yields simpler representations.

        std::vector< varPtr_t > normalisedVariants;
        for ( varPtr_t var : variants )
        {
            if ( var->isSNP() )
            {
                normalisedVariants.push_back( var );
                continue;
            }
            // Stores the current hit list of variants in reverse order. Will only contain one indel at any time.
            std::vector< varPtr_t > currentVars = {var};

            while ( not currentVars.empty() )
            {
                auto currentItem = currentVars.back();
                currentVars.pop_back();

                if ( currentItem->isSNP() )
                {
                    normalisedVariants.push_back( currentItem );
                }
                else if ( normalisedVariants.empty() )
                {
                    normalisedVariants.push_back( currentItem->getLeftAligned( referenceSequence.start() ) );
                }
                else
                {
                    currentItem = currentItem->getLeftAligned( normalisedVariants.back()->end() );

                    // Might have aligned next a variant we have seen before otherwise variant good to go.
                    if ( normalisedVariants.back()->joinable( currentItem ) )
                    {
                        // If so join, split and keep these to add to what needs to be normalized.
                        const auto newFrontElements = normalisedVariants.back()->join( currentItem )->split();
                        currentVars.insert( currentVars.end(), newFrontElements.rbegin(), newFrontElements.rend() );
                        normalisedVariants.pop_back();
                    }
                    else
                    {
                        normalisedVariants.push_back( currentItem );
                    }
                }
            }
        }
        return normalisedVariants;
    }

    //-----------------------------------------------------------------------------------------

    VariantGenerator::VariantGenerator( const utils::referenceSequencePtr_t & refSeq,
                                        const phred_t minBaseQual,
                                        const phred_t minMappingQual )
        : m_referenceSequence( refSeq ), m_minBaseQual( minBaseQual ), m_minMappingQual( minMappingQual )
    {
    }

    VariantContainer VariantGenerator::generateVariantsFromReads( io::perSampleRegionsReads_t perSamReadRanges ) const
    {
        VariantContainer variantContainer( m_minBaseQual, m_minMappingQual );

        for ( const auto & sampleAndReadRange : perSamReadRanges )
        {
            const auto & sampleName = sampleAndReadRange.first;
            for ( auto itRead = sampleAndReadRange.second.begin(); itRead != sampleAndReadRange.second.end(); ++itRead )
            {
                auto readPtr = itRead.getSharedPtr();

                const auto readReference = m_referenceSequence->subseq( readPtr->getRegion() );

                const auto readVariants = readPtr->getVariants();
                const auto normalisedVariants = normaliseVariantsOnStrand( readVariants, readReference );
                const auto readBreakpoints = readPtr->getBreakpoints();

                variantContainer.addVariantsFromRead( readPtr, normalisedVariants, readBreakpoints, sampleName );
            }
        }
        return variantContainer;
    }

    //-----------------------------------------------------------------------------------------
}
}
