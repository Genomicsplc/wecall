// All content Copyright (C) 2018 Genomics plc
#include "io/readfilters/readFilterAndTrimmer.hpp"
#include "utils/referenceSequence.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "caller/params.hpp"
#include "io/read.hpp"
using Cigar = wecall::alignment::Cigar;

BOOST_AUTO_TEST_CASE( testEmptyFilters )
{
    wecall::caller::params::Filters filterParams;
    filterParams.m_readMappingFilterQ = 0;
    filterParams.m_baseCallFilterN = 0;
    filterParams.m_baseCallFilterQ = 0;
    filterParams.m_duplicatesFilter = false;
    filterParams.m_noMatesFilter = false;
    filterParams.m_overlapTrim = false;
    filterParams.m_shortReadFilter = false;
    filterParams.m_shortReadTrim = false;

    wecall::io::ReadFilterAndTrimmer rft( filterParams );

    auto refSequence =
        std::make_shared< wecall::utils::ReferenceSequence >( wecall::caller::Region( "1", 100, 110 ), "ACGTAAAAGT" );
    wecall::io::Read::ReadParams readParams;
    readParams.m_sequenceStr = "WHY";
    readParams.m_qualities = "NOT";
    readParams.m_qname = "test";
    readParams.m_cigar = Cigar( "1M2I" );
    readParams.m_startPos = 100;
    readParams.m_endPos = 103;
    readParams.m_flag = BAM_FPROPER_PAIR;
    readParams.m_mappingQuality = 100;
    readParams.m_insertSize = 4;
    readParams.m_mateStartPos = 109;
    readParams.m_mateTid = 200;

    // test that a read with length passes trimAndFilter with all filters turned off
    auto emptyIntvlRead = std::make_shared< wecall::io::Read >( readParams, refSequence );
    BOOST_CHECK( rft.trimAndFilter( emptyIntvlRead ) );

    // test that a read with no length fails trimAndFilter
    readParams.m_cigar = Cigar( "3I" );
    readParams.m_startPos = 101;
    auto emptyIntvlReadZeroAlignedLength = std::make_shared< wecall::io::Read >( readParams, refSequence );
    BOOST_CHECK( not rft.trimAndFilter( emptyIntvlReadZeroAlignedLength ) );
}

BOOST_AUTO_TEST_CASE( testFilterFail )
{
    // define read which fails one of the hardcoded filters (isProperPair)
    wecall::caller::params::Filters filterParams;
    filterParams.m_readMappingFilterQ = 0;
    filterParams.m_baseCallFilterN = 0;
    filterParams.m_baseCallFilterQ = 0;
    filterParams.m_duplicatesFilter = false;
    filterParams.m_noMatesFilter = false;
    filterParams.m_overlapTrim = false;
    filterParams.m_shortReadFilter = false;
    filterParams.m_shortReadTrim = false;
    filterParams.m_allowImproperPairs = false;

    auto refSequence =
        std::make_shared< wecall::utils::ReferenceSequence >( wecall::caller::Region( "1", 100, 110 ), "ACGTAAAAGT" );

    wecall::io::ReadFilterAndTrimmer rft( filterParams );
    auto notProperPairRead = std::make_shared< wecall::io::Read >( "WHY", "NOT", "test", Cigar( "3M" ), 0, 101, 0, 100,
                                                                    200, 200, 0, refSequence );
    BOOST_CHECK( not notProperPairRead->isProperPair() );
    BOOST_CHECK( not rft.trimAndFilter( notProperPairRead ) );

    // add a filter
    filterParams.m_duplicatesFilter = true;  // filter all duplicates
    wecall::io::ReadFilterAndTrimmer rft_duplicates( filterParams );
    auto passesDuplicateFilterRead = std::make_shared< wecall::io::Read >(
        "WHY", "NOT", "test", Cigar( "3M" ), 0, 101, BAM_FPROPER_PAIR, 100, 200, 200, 0, refSequence );
    BOOST_CHECK( rft_duplicates.trimAndFilter( passesDuplicateFilterRead ) );
    auto failsDuplicateFilterRead = std::make_shared< wecall::io::Read >(
        "WHY", "NOT", "test", Cigar( "3M" ), 0, 101, BAM_FPROPER_PAIR + BAM_FDUP, 100, 200, 200, 0, refSequence );
    BOOST_CHECK( failsDuplicateFilterRead->isDuplicate() );
    BOOST_CHECK( not rft_duplicates.trimAndFilter( failsDuplicateFilterRead ) );
}

BOOST_AUTO_TEST_CASE( testOverlapTrimming )
{
    // Set overlap to last 1 base.  Qualities of the last base should be set to 0
    size_t startPos = 0;
    wecall::caller::params::Filters filterParams;
    filterParams.m_readMappingFilterQ = 0;
    filterParams.m_baseCallFilterN = 0;
    filterParams.m_baseCallFilterQ = 0;
    filterParams.m_duplicatesFilter = false;
    filterParams.m_noMatesFilter = false;
    filterParams.m_overlapTrim = true;
    filterParams.m_shortReadFilter = false;
    filterParams.m_shortReadTrim = false;

    std::string refSequenceStr = "AACGTAAAAGTGA";
    auto refSequence = std::make_shared< wecall::utils::ReferenceSequence >(
        wecall::caller::Region( "1", -3, -3 + refSequenceStr.size() ), refSequenceStr );

    wecall::io::ReadFilterAndTrimmer rftTrimOverlap( filterParams );
    auto overlappingRead = std::make_shared< wecall::io::Read >( "WHYWHYWHYA",  // seq
                                                                  "NOTNOTNOTB",  // qual
                                                                  "test",        // readGroupId
                                                                  Cigar( "10M" ), 0, startPos,
                                                                  BAM_FPROPER_PAIR + BAM_FMREVERSE,  // flag
                                                                  100,                               // mappingQual
                                                                  19,                                // insertSize
                                                                  8,                                 // mateStartPos
                                                                  0, refSequence );

    BOOST_CHECK( rftTrimOverlap.trimAndFilter( overlappingRead ) );
    BOOST_CHECK_EQUAL( overlappingRead->getQualities().substr( 0, 9 ), "NOTNOTNOT" );
    BOOST_CHECK_EQUAL( overlappingRead->getQualities().at( 9 ), constants::minAllowedQualityScore );
}

BOOST_AUTO_TEST_CASE( testShortReadTrimming )
{
    size_t startPos = 0;
    wecall::caller::params::Filters filterParams;
    filterParams.m_readMappingFilterQ = 0;
    filterParams.m_baseCallFilterN = 0;
    filterParams.m_baseCallFilterQ = 0;
    filterParams.m_noSimilarReadsFilter = false;
    filterParams.m_duplicatesFilter = false;
    filterParams.m_noMatesFilter = false;
    filterParams.m_overlapTrim = false;
    filterParams.m_shortReadFilter = false;
    filterParams.m_shortReadTrim = true;

    wecall::io::ReadFilterAndTrimmer rftTrimShortRead( filterParams );

    auto refSequenceStr = "ACGTAAAAGT" + std::string( 200, 'A' );
    auto refSequence = std::make_shared< wecall::utils::ReferenceSequence >(
        wecall::caller::Region( "1", -3, -3 + refSequenceStr.size() ), refSequenceStr );

    // test reverse short read
    auto reverseShortRead = std::make_shared< wecall::io::Read >( "WHYWHYWHYA",  // seq
                                                                   "NOTNOTNOTB",  // qual
                                                                   "test",        // readGroupId
                                                                   Cigar( "10M" ), 0, startPos,
                                                                   BAM_FPROPER_PAIR + BAM_FMREVERSE,  // flag
                                                                   100,                               // mappingQual
                                                                   9,                                 // insertSize
                                                                   8,                                 // mateStartPos
                                                                   0, refSequence );

    BOOST_CHECK( rftTrimShortRead.trimAndFilter( reverseShortRead ) );
    BOOST_CHECK_EQUAL( reverseShortRead->getQualities().substr( 0, 9 ), "NOTNOTNOT" );
    BOOST_CHECK_EQUAL( reverseShortRead->getQualities().at( 9 ), constants::minAllowedQualityScore );

    // test not reverse short read
    auto noReverseShortRead = std::make_shared< wecall::io::Read >( "WHYWHYWHYA",  // seq
                                                                     "NOTNOTNOTB",  // qual
                                                                     "test",        // readGroupId
                                                                     Cigar( "10M" ), 0, startPos,
                                                                     BAM_FPROPER_PAIR,  // flag
                                                                     100,               // mappingQual
                                                                     8,                 // insertSize
                                                                     8,                 // mateStartPos
                                                                     0, refSequence );

    BOOST_CHECK( rftTrimShortRead.trimAndFilter( noReverseShortRead ) );
    BOOST_CHECK_EQUAL( noReverseShortRead->getQualities().substr( 0, 8 ), "NOTNOTNO" );
    BOOST_CHECK_EQUAL( noReverseShortRead->getQualities().at( 8 ), constants::minAllowedQualityScore );
    BOOST_CHECK_EQUAL( noReverseShortRead->getQualities().at( 9 ), constants::minAllowedQualityScore );
}
