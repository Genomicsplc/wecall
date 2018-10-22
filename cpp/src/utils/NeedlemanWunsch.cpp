// All content Copyright (C) 2018 Genomics plc
#include "utils/NeedlemanWunsch.hpp"
//#include <boost/numeric/ublas/io.hpp>
#include <stack>
#include "utils/interval.hpp"

namespace wecall
{
namespace utils
{
    std::vector< NWVariant > NeedlemanWunsch::traceBack() const
    {
        std::stack< Backtrace > currentItems;

        const auto lastRefIndex = m_scoreMatrix.size1() - 1;
        const auto lastAltIndex = m_scoreMatrix.size2() - 1;
        currentItems.push(
            Backtrace{{}, lastRefIndex, lastAltIndex, m_traceMatrix( lastRefIndex, lastAltIndex ), false} );

        std::vector< std::vector< NWVariant > > allVariantCombinations;

        while ( not currentItems.empty() )
        {
            bool deadEnd = false;

            auto current = currentItems.top();
            currentItems.pop();

            while ( ( current.refIndex > 0 or current.altIndex > 0 ) and not deadEnd )
            {
                auto cache = current;
                if ( current.trace == 0 )
                {
                    // Reached the left side or the top
                    if ( current.refIndex and current.previousIsDeletion )
                    {
                        deadEnd = true;
                    }
                    else
                    {
                        if ( current.altIndex != 0 )
                        {
                            // insertion of the rest of the alt sequence
                            const auto alt = m_altString.substr( 0, current.altIndex );
                            current.variants.push_back( NWVariant( 0, 0, alt ) );
                        }
                        if ( current.refIndex != 0 )
                        {
                            // deletion of the rest of ref sequence
                            current.variants.push_back( NWVariant( 0, current.refIndex, "" ) );
                        }
                    }
                    break;
                }
                else if ( current.trace & m_matchBit )
                {
                    // No indel
                    current.trace ^= m_matchBit;
                    current.previousIsDeletion = false;
                    current.altIndex -= 1;
                    current.refIndex -= 1;
                }
                else if ( current.trace & m_snpBit )
                {
                    // SNP
                    current.trace ^= m_snpBit;
                    current.previousIsDeletion = false;
                    current.altIndex -= 1;
                    current.refIndex -= 1;
                    const auto alt = m_altString.substr( current.altIndex, 1 );
                    current.variants.push_back( NWVariant( current.refIndex, current.refIndex + 1, alt ) );
                }
                else if ( current.trace & m_insertSingleBit )
                {
                    // insertion of size 1
                    current.trace ^= m_insertSingleBit;
                    if ( current.previousIsDeletion )
                    {
                        deadEnd = true;
                    }
                    else
                    {
                        current.previousIsDeletion = false;
                        current.altIndex -= 1;
                        const auto alt = m_altString.substr( current.altIndex, 1 );
                        current.variants.push_back( NWVariant( current.refIndex, current.refIndex, alt ) );
                    }
                }
                else if ( current.trace & m_deleteSingleBit )
                {
                    // deletion of size 1
                    current.trace ^= m_deleteSingleBit;
                    current.previousIsDeletion = true;
                    current.refIndex -= 1;
                    current.variants.push_back( NWVariant( current.refIndex, current.refIndex + 1, "" ) );
                }
                else if ( current.trace & m_insertMultiBit )
                {
                    // insertion of size > 1
                    current.trace ^= m_insertMultiBit;
                    if ( current.previousIsDeletion )
                    {
                        deadEnd = true;
                    }
                    else
                    {
                        current.previousIsDeletion = false;
                        auto targetScore = m_scoreMatrix( current.refIndex, current.altIndex );
                        auto startingAltIndex = current.altIndex;

                        for ( std::size_t len = 2; len <= startingAltIndex; ++len )
                        {
                            current.altIndex = startingAltIndex - len;
                            auto actualScore = m_scoreMatrix( current.refIndex, current.altIndex ) +
                                               m_nWPenalties.insertionFunction( current.altIndex, len );

                            if ( targetScore == actualScore )
                            {
                                const auto alt = m_altString.substr( current.altIndex, len );
                                const auto variant = NWVariant( current.refIndex, current.refIndex, alt );
                                const auto newTrace = m_traceMatrix( current.refIndex, current.altIndex );

                                if ( newTrace == 0 )
                                {
                                    // We are at the beginning of the alt sequence, so there can't be any more starting
                                    // points that achieve the required score after this. Add the new variant to the
                                    // current variants and continue.
                                    current.variants.push_back( variant );
                                    break;
                                }
                                else
                                {
                                    // There may be more starting points earlier in the alt sequence, so we save the
                                    // information about starting at this point on the stack rather than continuing with
                                    // it straight away.
                                    auto newVariants = current.variants;
                                    newVariants.push_back( variant );
                                    currentItems.push( Backtrace{newVariants,
                                                                 current.refIndex,
                                                                 current.altIndex,
                                                                 newTrace,
                                                                 current.previousIsDeletion} );
                                }
                            }
                            if ( m_traceMatrix( current.refIndex, current.altIndex ) == 0 )
                            {
                                // If we get here then we got to the beginning of the alt sequence, and it wasn't a
                                // possible starting point for the insertion (if it was then we would have broken out of
                                // the for loop), so we shouldn't continue from the current position. We make this
                                // happen by setting deadEnd to true.
                                deadEnd = true;
                            }
                        }
                    }
                }
                else if ( current.trace & m_deleteMultiBit )
                {
                    // deletion of size > 1
                    current.trace ^= m_deleteMultiBit;
                    current.previousIsDeletion = true;

                    auto targetScore = m_scoreMatrix( current.refIndex, current.altIndex );
                    auto startingRefIndex = current.refIndex;

                    for ( std::size_t len = 2; len <= startingRefIndex; ++len )
                    {
                        current.refIndex = startingRefIndex - len;
                        auto actualScore = m_scoreMatrix( current.refIndex, current.altIndex ) +
                                           m_nWPenalties.deletionFunction( current.refIndex, len );

                        if ( targetScore == actualScore )
                        {
                            const auto variant = NWVariant( current.refIndex, startingRefIndex, "" );
                            const auto newTrace = m_traceMatrix( current.refIndex, current.altIndex );

                            if ( newTrace == 0 )
                            {
                                // We are at the beginning of the ref sequence, so there can't be any more starting
                                // points that achieve the required score after this. Add the new variant to the
                                // current variants and continue.
                                current.variants.push_back( variant );
                                break;
                            }
                            else
                            {
                                // There may be more starting points earlier in the ref sequence, so we save the
                                // information about starting at this point on the stack rather than continuing with
                                // it straight away.
                                auto newVariants = current.variants;
                                newVariants.push_back( variant );
                                currentItems.push( Backtrace{newVariants,
                                                             current.refIndex,
                                                             current.altIndex,
                                                             newTrace,
                                                             current.previousIsDeletion} );
                            }
                        }
                        if ( m_traceMatrix( current.refIndex, current.altIndex ) == 0 )
                        {
                            // If we get here then we got to the beginning of the ref sequence, and it wasn't a
                            // possible starting point for the insertion (if it was then we would have broken out of
                            // the for loop), so we shouldn't continue from the current position. We make this
                            // happen by setting deadEnd to true.
                            deadEnd = true;
                        }
                    }
                }
                if ( current.trace != 0 )
                {
                    cache.trace = current.trace;
                    currentItems.push( cache );
                }
                current.trace = m_traceMatrix( current.refIndex, current.altIndex );
            }
            if ( not deadEnd )
            {
                allVariantCombinations.push_back( current.variants );
            }
        }
        if ( allVariantCombinations.size() > 1 )
        {
            ECHIDNA_LOG( DEBUG, "Warning: while normalizing variants we found "
                                    << allVariantCombinations.size() << " possible normalizations" << std::endl );
        }

        return allVariantCombinations.back();
    }

    int32_t NeedlemanWunsch::getScoreMatrix()
    {
        const std::size_t refDim = m_referenceString.size() + 1;
        const std::size_t altDim = m_altString.size() + 1;

        //        scoreMatrix_t m_scoreMatrix( refDim, altDim );
        //        traceMatrix_t m_traceMatrix( refDim, altDim );
        m_scoreMatrix = scoreMatrix_t( refDim, altDim );
        m_traceMatrix = traceMatrix_t( refDim, altDim );

        m_scoreMatrix( 0, 0 ) = 0;

        // Initialise the first column of the matrix
        //        for ( std::size_t rowIndex = 1; rowIndex < refDim; ++rowIndex )
        //        {
        //            m_scoreMatrix( rowIndex, 0 ) = m_nWPenalties.insertionFunction(0, rowIndex);
        //            m_traceMatrix( rowIndex, 0 ) = 0;
        //        }

        if ( refDim > 1 )
        {
            m_scoreMatrix( 1, 0 ) = m_nWPenalties.insertionOpen();
            m_traceMatrix( 1, 0 ) = 0;
        }

        for ( std::size_t rowIndex = 2; rowIndex < refDim; ++rowIndex )
        {
            m_scoreMatrix( rowIndex, 0 ) = m_scoreMatrix( rowIndex - 1, 0 ) + m_nWPenalties.extendInsertion();
            m_traceMatrix( rowIndex, 0 ) = 0;
        }

        // Initialise the first row of the matrix
        //        for ( std::size_t colIndex = 1; colIndex < refDim; ++colIndex )
        //        {
        //            m_scoreMatrix( 0, colIndex ) = m_nWPenalties.deletionFunction(0, colIndex);
        //            m_traceMatrix( 0, colIndex ) = 0;
        //        }

        if ( altDim > 1 )
        {
            m_scoreMatrix( 0, 1 ) = m_nWPenalties.deletionOpen();
            m_traceMatrix( 0, 1 ) = 0;
        }

        for ( std::size_t colIndex = 2; colIndex < altDim; ++colIndex )
        {
            m_scoreMatrix( 0, colIndex ) = m_scoreMatrix( 0, colIndex - 1 ) + m_nWPenalties.extendDeletion();
            m_traceMatrix( 0, colIndex ) = 0;
        }

        // Create the colScore matrix we will use when filling in the matrix
        // It took me a while to work out why we use these values. When we add
        std::vector< int32_t > colScore( altDim );
        // Note that colScore[0] is never used, so there's no need to initialise it
        //        for ( std::size_t colIndex = 1; colIndex < altDim; ++colIndex )
        //        {
        //            // When we add n*m_nWPenalties.gapLinear to colScore[colIndex] we want to get the score for having
        //            // an insertion of size colIndex followed by a deletion of size n
        //            colScore[colIndex] = m_nWPenalties.insertionFunction(0, colIndex) +
        //                                 m_nWPenalties.deletionFunction(0, 1) - m_nWPenalties.gapLinear;
        //        }
        if ( altDim > 1 )
        {
            colScore[1] = m_nWPenalties.insertionOpen() + m_nWPenalties.deletionOpen() - m_nWPenalties.extendDeletion();
        }

        for ( std::size_t colIndex = 2; colIndex < altDim; ++colIndex )
        {
            colScore[colIndex] = colScore[colIndex - 1] + m_nWPenalties.extendDeletion();
        }

        // Now go through the rest of the matrix, row by row, calculating the correct values
        for ( int32_t rowIndex = 1; rowIndex < static_cast< int32_t >( refDim ); ++rowIndex )
        {
            //            // When we add m_nWPenalties.gapLinear to rowScore we want to get the score for having a
            //            // deletion of size rowIndex followed by an insertion of size 1
            auto rowScore = m_nWPenalties.deletionFunction( 0, rowIndex ) + m_nWPenalties.insertionFunction( 0, 1 ) -
                            m_nWPenalties.extendInsertion();
            //            auto rowScore = 2 * m_nWPenalties.gapConstant + rowIndex * m_nWPenalties.extendDeletion();

            for ( int32_t colIndex = 1; colIndex < static_cast< int32_t >( altDim ); ++colIndex )
            {
                // Calculate score for not having a gap, i.e. sequences match or there is a SNP
                const auto refBase = m_referenceString.at( rowIndex - 1 );
                const auto altBase = m_altString.at( colIndex - 1 );
                const int32_t noGapScore =
                    m_scoreMatrix( rowIndex - 1, colIndex - 1 ) + m_nWPenalties.matchFunction( refBase, altBase );

                const auto insertSingleScore =
                    m_scoreMatrix( rowIndex, colIndex - 1 ) + m_nWPenalties.insertionFunction( colIndex - 1, 1 );
                const auto insertMultiScore = rowScore + m_nWPenalties.extendInsertion();
                rowScore = std::max( insertSingleScore, insertMultiScore );

                const auto deleteSingleScore =
                    m_scoreMatrix( rowIndex - 1, colIndex ) + m_nWPenalties.deletionFunction( rowIndex - 1, 1 );
                const auto deleteMultiScore = colScore[colIndex] + m_nWPenalties.extendDeletion();
                colScore[colIndex] = std::max( deleteSingleScore, deleteMultiScore );

                const auto bestScore =
                    std::max( {noGapScore, insertSingleScore, insertMultiScore, deleteSingleScore, deleteMultiScore} );

                m_scoreMatrix( rowIndex, colIndex ) = bestScore;

                bitType_t trace = 0;

                if ( bestScore == noGapScore )
                {
                    if ( refBase == altBase )
                    {
                        trace |= m_matchBit;
                    }
                    else
                    {
                        trace |= m_snpBit;
                    }
                }
                if ( bestScore == insertSingleScore )
                {
                    trace |= m_insertSingleBit;
                }
                if ( bestScore == deleteSingleScore )
                {
                    trace |= m_deleteSingleBit;
                }
                if ( bestScore == insertMultiScore )
                {
                    trace |= m_insertMultiBit;
                }
                if ( bestScore == deleteMultiScore )
                {
                    trace |= m_deleteMultiBit;
                }

                m_traceMatrix( rowIndex, colIndex ) = trace;
            }
        }

        return m_scoreMatrix( refDim - 1, altDim - 1 );
    }

    std::ostream & operator<<( std::ostream & out, const NWPenalties & nWPenalties )
    {
        out << "NWPenalties(" << nWPenalties.m_baseMatch << ", " << nWPenalties.m_baseMismatch << ", "
            << nWPenalties.m_insertionConstant << ", " << nWPenalties.m_insertionLinear << ", "
            << nWPenalties.m_deletionConstant << ", " << nWPenalties.m_deletionLinear << ", "
            << nWPenalties.m_indelPosition << ")";
        return out;
    }

    std::ostream & operator<<( std::ostream & out, const NWVariant & nWVariant )
    {
        out << "NWVariant(" << nWVariant.m_start << ", " << nWVariant.m_end << ", " << nWVariant.m_alt << ")";
        return out;
    }
}
}