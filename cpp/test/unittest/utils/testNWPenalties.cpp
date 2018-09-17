#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "utils/NeedlemanWunsch.hpp"

using echidna::utils::NWPenalties;

BOOST_AUTO_TEST_CASE( shouldMakeSameObjectWithBothInitialisers )
{
    const int32_t baseMatch = 1;
    const int32_t baseMismatch = 2;
    const int32_t indelConstant = 3;
    const int32_t indelLinear = 4;
    const int32_t indelPosition = 5;
    const auto nWPenalties1 =
        NWPenalties( baseMatch, baseMismatch, indelConstant, indelLinear, indelConstant, indelLinear, indelPosition );
    const auto nWPenalties2 = NWPenalties( baseMatch, baseMismatch, indelConstant, indelLinear, indelPosition );

    BOOST_CHECK_EQUAL( nWPenalties1, nWPenalties2 );
}

BOOST_AUTO_TEST_CASE( shouldCalculateMatchingBasesPenalties )
{
    const int32_t baseMatch = 1;
    const int32_t baseMismatch = 10;
    const int32_t insertionConstant = 10;
    const int32_t insertionLinear = 10;
    const int32_t deletionConstant = 10;
    const int32_t deletionLinear = 10;
    const int32_t indelPosition = 10;
    const auto nWPenalties = NWPenalties( baseMatch, baseMismatch, insertionConstant, insertionLinear, deletionConstant,
                                          deletionLinear, indelPosition );

    const char bases[4] = {'A', 'C', 'G', 'T'};

    for ( auto index = 0; index < 4; ++index )
    {
        BOOST_CHECK_EQUAL( baseMatch, nWPenalties.matchFunction( bases[index], bases[index] ) );
    }
}

BOOST_AUTO_TEST_CASE( shouldCalculateMismatchingBasesPenalties )
{
    const int32_t baseMatch = 10;
    const int32_t baseMismatch = 1;
    const int32_t insertionConstant = 10;
    const int32_t insertionLinear = 10;
    const int32_t deletionConstant = 10;
    const int32_t deletionLinear = 10;
    const int32_t indelPosition = 10;
    const auto nWPenalties = NWPenalties( baseMatch, baseMismatch, insertionConstant, insertionLinear, deletionConstant,
                                          deletionLinear, indelPosition );

    const char bases[4] = {'A', 'C', 'G', 'T'};

    for ( auto index1 = 0; index1 < 4; ++index1 )
    {
        for ( auto index2 = 0; index2 < 4; ++index2 )
        {
            if ( index1 != index2 )
            {
                BOOST_CHECK_EQUAL( baseMismatch, nWPenalties.matchFunction( bases[index1], bases[index2] ) );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( shouldCalculateInsertionOpenPenalty )
{
    const int32_t baseMatch = 10;
    const int32_t baseMismatch = 10;
    const int32_t insertionConstant = 1;
    const int32_t insertionLinear = 2;
    const int32_t deletionConstant = 10;
    const int32_t deletionLinear = 10;
    const int32_t indelPosition = 10;
    const auto nWPenalties = NWPenalties( baseMatch, baseMismatch, insertionConstant, insertionLinear, deletionConstant,
                                          deletionLinear, indelPosition );

    BOOST_CHECK_EQUAL( insertionConstant + insertionLinear, nWPenalties.insertionOpen() );
}

BOOST_AUTO_TEST_CASE( shouldCalculateDeletionOpenPenalty )
{
    const int32_t baseMatch = 10;
    const int32_t baseMismatch = 10;
    const int32_t insertionConstant = 10;
    const int32_t insertionLinear = 10;
    const int32_t deletionConstant = 1;
    const int32_t deletionLinear = 2;
    const int32_t indelPosition = 10;
    const auto nWPenalties = NWPenalties( baseMatch, baseMismatch, insertionConstant, insertionLinear, deletionConstant,
                                          deletionLinear, indelPosition );

    BOOST_CHECK_EQUAL( deletionConstant + deletionLinear, nWPenalties.deletionOpen() );
}

BOOST_AUTO_TEST_CASE( shouldCalculateInsertionFunction )
{
    const int32_t baseMatch = 10;
    const int32_t baseMismatch = 10;
    const int32_t insertionConstant = 1;
    const int32_t insertionLinear = 2;
    const int32_t deletionConstant = 10;
    const int32_t deletionLinear = 10;
    const int32_t indelPosition = 10;
    const auto nWPenalties = NWPenalties( baseMatch, baseMismatch, insertionConstant, insertionLinear, deletionConstant,
                                          deletionLinear, indelPosition );

    for ( auto position = 0; position < 10; ++position )
    {
        BOOST_CHECK_EQUAL( 0, nWPenalties.insertionFunction( position, 0 ) );

        for ( auto insertionLength = 1; insertionLength < 10; ++insertionLength )
        {
            BOOST_CHECK_EQUAL( insertionConstant + insertionLength * insertionLinear + position * indelPosition,
                               nWPenalties.insertionFunction( position, insertionLength ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( shouldCalculateDeletionFunction )
{
    const int32_t baseMatch = 10;
    const int32_t baseMismatch = 10;
    const int32_t insertionConstant = 10;
    const int32_t insertionLinear = 10;
    const int32_t deletionConstant = 1;
    const int32_t deletionLinear = 2;
    const int32_t indelPosition = 10;
    const auto nWPenalties = NWPenalties( baseMatch, baseMismatch, insertionConstant, insertionLinear, deletionConstant,
                                          deletionLinear, indelPosition );

    for ( auto position = 0; position < 10; ++position )
    {
        BOOST_CHECK_EQUAL( 0, nWPenalties.insertionFunction( position, 0 ) );

        for ( auto deletionLength = 1; deletionLength < 10; ++deletionLength )
        {
            BOOST_CHECK_EQUAL( deletionConstant + deletionLength * deletionLinear + position * indelPosition,
                               nWPenalties.deletionFunction( position, deletionLength ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( shouldCalculateExtendInsertion )
{
    const int32_t baseMatch = 10;
    const int32_t baseMismatch = 10;
    const int32_t insertionConstant = 1;
    const int32_t insertionLinear = 2;
    const int32_t deletionConstant = 10;
    const int32_t deletionLinear = 10;
    const int32_t indelPosition = 10;
    const auto nWPenalties = NWPenalties( baseMatch, baseMismatch, insertionConstant, insertionLinear, deletionConstant,
                                          deletionLinear, indelPosition );

    for ( auto position = 0; position < 10; ++position )
    {
        for ( auto insertionLength = 1; insertionLength < 10; ++insertionLength )
        {
            BOOST_CHECK_EQUAL(
                nWPenalties.insertionFunction( position, insertionLength + 1 ),
                nWPenalties.insertionFunction( position, insertionLength ) + nWPenalties.extendInsertion() );
        }
    }
}

BOOST_AUTO_TEST_CASE( shouldCalculateExtendDeletion )
{
    const int32_t baseMatch = 10;
    const int32_t baseMismatch = 10;
    const int32_t insertionConstant = 10;
    const int32_t insertionLinear = 10;
    const int32_t deletionConstant = 1;
    const int32_t deletionLinear = 2;
    const int32_t indelPosition = 10;
    const auto nWPenalties = NWPenalties( baseMatch, baseMismatch, insertionConstant, insertionLinear, deletionConstant,
                                          deletionLinear, indelPosition );

    for ( auto position = 0; position < 10; ++position )
    {
        for ( auto deletionLength = 1; deletionLength < 10; ++deletionLength )
        {
            BOOST_CHECK_EQUAL(
                nWPenalties.deletionFunction( position, deletionLength + 1 ),
                nWPenalties.deletionFunction( position, deletionLength ) + nWPenalties.extendDeletion() );
        }
    }
}