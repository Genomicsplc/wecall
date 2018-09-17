#include "alignment/cigar.hpp"
#include "samtools/bam.h"
#include "alignment/cigarItems.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <vector>

using echidna::alignment::Cigar;
using echidna::alignment::CigarItem;
using echidna::alignment::cigarFlags;
using echidna::alignment::CigarItemBuilder;

BOOST_AUTO_TEST_CASE( shouldRollFrom0To6 )
{
    std::vector< std::string > cigarStrings;
    for ( uint32_t cigarInt = 0; cigarInt < 9u; ++cigarInt )
    {
        cigarStrings.emplace_back( CigarItemBuilder::roll( cigarInt )->toString() );
    }
    std::vector< std::string > expectedCigarItems = {"0M", "0I", "0D", "0N", "0S", "0H", "0P", "0=", "0X"};
    BOOST_CHECK_EQUAL_COLLECTIONS( cigarStrings.begin(), cigarStrings.end(), expectedCigarItems.begin(),
                                   expectedCigarItems.end() );
}

BOOST_AUTO_TEST_CASE( shouldRollMatch )
{
    auto cigar = CigarItemBuilder::roll( cigarFlags::MATCH, 2 );
    BOOST_CHECK( dynamic_cast< echidna::alignment::CigarMatch * >( cigar.get() ) != nullptr );
}

BOOST_AUTO_TEST_CASE( shouldRollDeletion )
{
    auto cigar = CigarItemBuilder::roll( cigarFlags::DELETION, 2 );
    BOOST_CHECK( dynamic_cast< echidna::alignment::CigarDeletion * >( cigar.get() ) != nullptr );
}

BOOST_AUTO_TEST_CASE( shouldRollInsertion )
{
    auto cigar = CigarItemBuilder::roll( cigarFlags::INSERTION, 2 );
    BOOST_CHECK( dynamic_cast< echidna::alignment::CigarInsertion * >( cigar.get() ) != nullptr );
}

BOOST_AUTO_TEST_CASE( shouldRollSkip )
{
    auto cigar = CigarItemBuilder::roll( cigarFlags::SKIP, 2 );
    BOOST_CHECK( dynamic_cast< echidna::alignment::CigarSkip * >( cigar.get() ) != nullptr );
}

BOOST_AUTO_TEST_CASE( shouldRollHardClip )
{
    auto cigar = CigarItemBuilder::roll( cigarFlags::HARD_CLIP, 2 );
    BOOST_CHECK( dynamic_cast< echidna::alignment::CigarHardClip * >( cigar.get() ) != nullptr );
}

BOOST_AUTO_TEST_CASE( shouldRollSoftClip )
{
    auto cigar = CigarItemBuilder::roll( cigarFlags::SOFT_CLIP, 2 );
    BOOST_CHECK( dynamic_cast< echidna::alignment::CigarSoftClip * >( cigar.get() ) != nullptr );
}

BOOST_AUTO_TEST_CASE( shouldRollPadding )
{
    auto cigar = CigarItemBuilder::roll( cigarFlags::PADDING, 2 );
    BOOST_CHECK( dynamic_cast< echidna::alignment::CigarPadding * >( cigar.get() ) != nullptr );
}

BOOST_AUTO_TEST_CASE( shouldComputeLengthsOneMatchItem )
{
    Cigar cigar( "4M" );
    BOOST_CHECK_EQUAL( cigar.length(), 4 );
    BOOST_CHECK_EQUAL( cigar.lengthInRef(), 4 );
}

BOOST_AUTO_TEST_CASE( shouldComputeLengthsMultipleItems )
{
    Cigar cigar( "4M7I" );
    BOOST_CHECK_EQUAL( cigar.length(), 11 );
    BOOST_CHECK_EQUAL( cigar.lengthInRef(), 4 );
}

BOOST_AUTO_TEST_CASE( shouldComputeLengthWithoutSoftCippling )
{
    Cigar cigar( "2S4M1S4S7I100S" );
    BOOST_CHECK_EQUAL( cigar.lengthInSeqWithoutSoftClipping(), 11 );
}

BOOST_AUTO_TEST_CASE( shouldGetTheCorrectRelativePositions )
{
    using echidna::alignment::emptyPos;
    Cigar cigar( "2M1I3D1M" );
    auto actualPositions = cigar.getRefPositions( 0 );
    std::vector< int64_t > expectedPositions{0, 1, emptyPos, 5};
    BOOST_CHECK_EQUAL_COLLECTIONS( actualPositions.begin(), actualPositions.end(), expectedPositions.begin(),
                                   expectedPositions.end() );
}

BOOST_AUTO_TEST_CASE( testShouldOutputCorrectCigarString )
{
    std::string cigarString = "2M1I3D4P5H8S9N123567=342342X";
    Cigar cigar( cigarString );
    BOOST_CHECK_EQUAL( cigarString, cigar.toString() );
}

BOOST_AUTO_TEST_CASE( shouldGetEmptyIntervalIfInputIsEmptyAndCigarMatch )
{
    Cigar cigar( "5M" );
    echidna::utils::Interval inputInterval( 1L, 1L );
    echidna::utils::Interval expectedResult( 1L, 1L );
    BOOST_CHECK_EQUAL( expectedResult, cigar.getInverseInterval( inputInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldGetMatchingIntervalIfInputIsNonEmptyAndCigarMatch )
{
    Cigar cigar( "5M" );
    echidna::utils::Interval inputInterval( 1L, 4L );
    echidna::utils::Interval expectedResult( 1L, 4L );
    BOOST_CHECK_EQUAL( expectedResult, cigar.getInverseInterval( inputInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldGetIntervalCorrespondingToInsertion )
{
    Cigar cigar( "1M2I3M" );
    echidna::utils::Interval inputInterval( 1L, 1L );
    echidna::utils::Interval expectedResult( 1L, 3L );
    BOOST_CHECK_EQUAL( expectedResult, cigar.getInverseInterval( inputInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldGetIntervalCorrespondingToDeletion )
{
    Cigar cigar( "1M2D3M" );
    echidna::utils::Interval inputInterval( 1L, 3L );
    echidna::utils::Interval expectedResult( 1L, 1L );
    BOOST_CHECK_EQUAL( expectedResult, cigar.getInverseInterval( inputInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldGetIntervalCorrespondingToMatchAfterInsertion )
{
    Cigar cigar( "1M2I3M" );
    echidna::utils::Interval inputInterval( 2L, 3L );
    echidna::utils::Interval expectedResult( 4L, 5L );
    BOOST_CHECK_EQUAL( expectedResult, cigar.getInverseInterval( inputInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldGetIntervalCorrespondingToMatchAfterDeletion )
{
    Cigar cigar( "1M2D3M" );
    echidna::utils::Interval inputInterval( 4L, 5L );
    echidna::utils::Interval expectedResult( 2L, 3L );
    BOOST_CHECK_EQUAL( expectedResult, cigar.getInverseInterval( inputInterval ) );
}

BOOST_AUTO_TEST_CASE( shouldClipCigarWithOnlySoftClipping )
{
    Cigar cigar( "10S" );

    const auto clippings = cigar.stripSoftClipping();
    BOOST_CHECK_EQUAL( std::distance( cigar.cbegin(), cigar.cend() ), 0 );
    BOOST_CHECK_EQUAL( cigar.lengthInSeq(), 0L );
    BOOST_CHECK_EQUAL( cigar.lengthInRef(), 0L );
    BOOST_CHECK_EQUAL( clippings.first, 10L );
    BOOST_CHECK_EQUAL( clippings.second, 0L );
}

BOOST_AUTO_TEST_CASE( shouldClipEmptyCigar )
{
    Cigar cigar( "" );

    const auto clippings = cigar.stripSoftClipping();
    BOOST_CHECK_EQUAL( std::distance( cigar.cbegin(), cigar.cend() ), 0 );
    BOOST_CHECK_EQUAL( cigar.lengthInSeq(), 0L );
    BOOST_CHECK_EQUAL( cigar.lengthInRef(), 0L );
    BOOST_CHECK_EQUAL( clippings.first, 0L );
    BOOST_CHECK_EQUAL( clippings.second, 0L );
}

BOOST_AUTO_TEST_CASE( shouldClipSoftClippingAtStart )
{
    Cigar cigar( "10S2M" );

    const auto clippings = cigar.stripSoftClipping();
    BOOST_CHECK_EQUAL( std::distance( cigar.cbegin(), cigar.cend() ), 1 );
    BOOST_CHECK_EQUAL( cigar.lengthInSeq(), 2L );
    BOOST_CHECK_EQUAL( cigar.lengthInRef(), 2L );
    BOOST_CHECK_EQUAL( clippings.first, 10L );
    BOOST_CHECK_EQUAL( clippings.second, 0L );
}

BOOST_AUTO_TEST_CASE( shouldClipSoftClippingAtEnd )
{
    Cigar cigar( "2M10S" );

    const auto clippings = cigar.stripSoftClipping();
    BOOST_CHECK_EQUAL( std::distance( cigar.cbegin(), cigar.cend() ), 1 );
    BOOST_CHECK_EQUAL( cigar.lengthInSeq(), 2L );
    BOOST_CHECK_EQUAL( cigar.lengthInRef(), 2L );
    BOOST_CHECK_EQUAL( clippings.first, 0L );
    BOOST_CHECK_EQUAL( clippings.second, 10L );
}

BOOST_AUTO_TEST_CASE( shouldClipSoftClippingAtStartAndEnd )
{
    Cigar cigar( "8S2M10S" );

    const auto clippings = cigar.stripSoftClipping();
    BOOST_CHECK_EQUAL( std::distance( cigar.cbegin(), cigar.cend() ), 1 );
    BOOST_CHECK_EQUAL( cigar.lengthInSeq(), 2L );
    BOOST_CHECK_EQUAL( cigar.lengthInRef(), 2L );
    BOOST_CHECK_EQUAL( clippings.first, 8L );
    BOOST_CHECK_EQUAL( clippings.second, 10L );
}

BOOST_AUTO_TEST_CASE( shouldNotClipInMiddle )
{
    Cigar cigar( "2M10S2M" );

    const auto clippings = cigar.stripSoftClipping();
    BOOST_CHECK_EQUAL( std::distance( cigar.cbegin(), cigar.cend() ), 3 );

    BOOST_CHECK_EQUAL( clippings.first, 0L );
    BOOST_CHECK_EQUAL( clippings.second, 0L );
}
