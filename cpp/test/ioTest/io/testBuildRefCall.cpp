// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "io/read.hpp"
#include "io/readIntervalTree.hpp"
#include "io/readRange.hpp"
#include "io/readSummaries.hpp"
#include "caller/diploid/referenceCalling.hpp"
#include "vcf/field.hpp"
#include <cmath>

using Read = echidna::io::Read;
using Cigar = echidna::alignment::Cigar;
using echidna::caller::Region;
using echidna::caller::Call;
using echidna::io::RegionsReads;
using echidna::utils::BasePairSequence;
using echidna::caller::model::buildRefCall;
using Annotation = echidna::caller::Annotation;
using echidna::vcf::info::DP_key;

BOOST_AUTO_TEST_CASE( shouldCallRefFor1ReadAnd1Sample )
{
    Region region = Region( "1", 0, 5 );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( region, std::string( 5, 'A' ) );

    const int64_t startPos = 0;
    const auto read1 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    echidna::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );

    RegionsReads regionSetReads( region, readContainer.getFullRange(), 0 );

    echidna::io::perSampleRegionsReads_t perSampleReads = {{"sample1", regionSetReads}};

    const auto refCalls = buildRefCall( region, perSampleReads, 10, {2}, 0.2 );

    BOOST_REQUIRE_EQUAL( refCalls.size(), 2 );

    BOOST_CHECK_EQUAL( refCalls[0].getAnnotation( Annotation::BEG ), 1 );
    BOOST_CHECK_EQUAL( refCalls[0].getAnnotation( Annotation::END ), 4 );
    BOOST_CHECK_EQUAL( refCalls[0].getAnnotation( Annotation::LEN ), 4 );
    BOOST_CHECK_EQUAL( refCalls[0].samples[0].getAnnotation( Annotation::MIN_DP ), 1 );
    BOOST_CHECK_EQUAL( refCalls[0].samples[0].getAnnotation( Annotation::FORMAT_DP ), 1 );
    BOOST_CHECK_CLOSE( refCalls[0].samples[0].getAnnotation( Annotation::GQ ), 3.0103, 1 );

    BOOST_CHECK_EQUAL( refCalls[1].getAnnotation( Annotation::BEG ), 5 );
    BOOST_CHECK_EQUAL( refCalls[1].getAnnotation( Annotation::END ), 5 );
    BOOST_CHECK_EQUAL( refCalls[1].getAnnotation( Annotation::LEN ), 1 );
    BOOST_CHECK_EQUAL( refCalls[1].samples[0].getAnnotation( Annotation::MIN_DP ), 0 );
    BOOST_CHECK_EQUAL( refCalls[1].samples[0].getAnnotation( Annotation::FORMAT_DP ), 0 );
    BOOST_CHECK_EQUAL( std::isnan( refCalls[1].samples[0].getAnnotation( Annotation::GQ ) ), true );
}

BOOST_AUTO_TEST_CASE( shouldCallRefFor2OverlappingReadsAnd1Sample )
{
    Region region = Region( "1", 0, 5 );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( region, std::string( 5, 'A' ) );

    const int64_t startPos = 0;
    const auto read1 = std::make_shared< Read >( BasePairSequence( 3, 'A' ), std::string( 3, 'Q' ), "0", Cigar( "3M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const auto read2 = std::make_shared< Read >( BasePairSequence( 3, 'A' ), std::string( 3, 'Q' ), "0", Cigar( "3M" ),
                                                 0, startPos + 2, 0, 0, 0, 0, 0, refSequence );

    echidna::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );
    readContainer.insert( read2 );

    RegionsReads regionSetReads( region, readContainer.getFullRange(), 0 );

    echidna::io::perSampleRegionsReads_t perSampleReads = {{"sample1", regionSetReads}};

    const auto refCalls = buildRefCall( region, perSampleReads, 10, {2}, 0.2 );

    BOOST_REQUIRE_EQUAL( refCalls.size(), 3 );

    BOOST_CHECK_EQUAL( refCalls[0].getAnnotation( Annotation::BEG ), 1 );
    BOOST_CHECK_EQUAL( refCalls[0].getAnnotation( Annotation::END ), 2 );
    BOOST_CHECK_EQUAL( refCalls[0].getAnnotation( Annotation::LEN ), 2 );
    BOOST_CHECK_EQUAL( refCalls[0].samples[0].getAnnotation( Annotation::MIN_DP ), 1 );
    BOOST_CHECK_EQUAL( refCalls[0].samples[0].getAnnotation( Annotation::FORMAT_DP ), 1 );
    BOOST_CHECK_CLOSE( refCalls[0].samples[0].getAnnotation( Annotation::GQ ), 3.0103, 1 );

    BOOST_CHECK_EQUAL( refCalls[1].getAnnotation( Annotation::BEG ), 3 );
    BOOST_CHECK_EQUAL( refCalls[1].getAnnotation( Annotation::END ), 3 );
    BOOST_CHECK_EQUAL( refCalls[1].getAnnotation( Annotation::LEN ), 1 );
    BOOST_CHECK_EQUAL( refCalls[1].samples[0].getAnnotation( Annotation::MIN_DP ), 2 );
    BOOST_CHECK_EQUAL( refCalls[1].samples[0].getAnnotation( Annotation::FORMAT_DP ), 2 );
    BOOST_CHECK_CLOSE( refCalls[1].samples[0].getAnnotation( Annotation::GQ ), 5.9159, 1 );

    BOOST_CHECK_EQUAL( refCalls[2].getAnnotation( Annotation::BEG ), 4 );
    BOOST_CHECK_EQUAL( refCalls[2].getAnnotation( Annotation::END ), 5 );
    BOOST_CHECK_EQUAL( refCalls[2].getAnnotation( Annotation::LEN ), 2 );
    BOOST_CHECK_EQUAL( refCalls[2].samples[0].getAnnotation( Annotation::MIN_DP ), 1 );
    BOOST_CHECK_EQUAL( refCalls[2].samples[0].getAnnotation( Annotation::FORMAT_DP ), 1 );
    BOOST_CHECK_CLOSE( refCalls[2].samples[0].getAnnotation( Annotation::GQ ), 3.0103, 1 );
}

BOOST_AUTO_TEST_CASE( shouldCallRefForManyReadsFor1Sample )
{
    Region region = Region( "1", 0, 5 );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( region, std::string( 5, 'A' ) );

    const int64_t startPos = 0;
    const auto read1 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const auto read2 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const auto read3 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const auto read4 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const auto read5 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const auto read6 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const auto read7 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const auto read8 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const auto read9 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const auto read10 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                  0, startPos + 1, 0, 0, 0, 0, 0, refSequence );

    echidna::io::readIntervalTree_t readContainer( 0, 100 );
    readContainer.insert( read1 );
    readContainer.insert( read2 );
    readContainer.insert( read3 );
    readContainer.insert( read4 );
    readContainer.insert( read5 );
    readContainer.insert( read6 );
    readContainer.insert( read7 );
    readContainer.insert( read8 );
    readContainer.insert( read9 );
    readContainer.insert( read10 );

    RegionsReads regionSetReads( region, readContainer.getFullRange(), 0 );

    echidna::io::perSampleRegionsReads_t perSampleReads = {{"sample1", regionSetReads}};

    const auto refCalls = buildRefCall( region, perSampleReads, 10, {2}, 0.1 );

    BOOST_REQUIRE_EQUAL( refCalls.size(), 2 );

    BOOST_CHECK_EQUAL( refCalls[0].getAnnotation( Annotation::BEG ), 1 );
    BOOST_CHECK_EQUAL( refCalls[0].getAnnotation( Annotation::END ), 4 );
    BOOST_CHECK_EQUAL( refCalls[0].getAnnotation( Annotation::LEN ), 4 );
    BOOST_CHECK_EQUAL( refCalls[0].samples[0].getAnnotation( Annotation::MIN_DP ), 9 );
    BOOST_CHECK_EQUAL( refCalls[0].samples[0].getAnnotation( Annotation::FORMAT_DP ), 10 );  // rounded 9.75
    BOOST_CHECK_CLOSE( refCalls[0].samples[0].getAnnotation( Annotation::GQ ), 23.85277, 1 );

    BOOST_CHECK_EQUAL( refCalls[1].getAnnotation( Annotation::BEG ), 5 );
    BOOST_CHECK_EQUAL( refCalls[1].getAnnotation( Annotation::END ), 5 );
    BOOST_CHECK_EQUAL( refCalls[1].getAnnotation( Annotation::LEN ), 1 );
    BOOST_CHECK_EQUAL( refCalls[1].samples[0].getAnnotation( Annotation::MIN_DP ), 1 );
    BOOST_CHECK_EQUAL( refCalls[1].samples[0].getAnnotation( Annotation::FORMAT_DP ), 1 );
    BOOST_CHECK_CLOSE( refCalls[1].samples[0].getAnnotation( Annotation::GQ ), 3.0103, 1 );
}

BOOST_AUTO_TEST_CASE( shouldCallRefFor1ReadsFor2Samples )
{
    Region region = Region( "1", 0, 5 );
    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( region, std::string( 5, 'A' ) );

    const int64_t startPos = 0;
    const auto read1 = std::make_shared< Read >( BasePairSequence( 3, 'A' ), std::string( 3, 'Q' ), "0", Cigar( "3M" ),
                                                 0, startPos, 0, 0, 0, 0, 0, refSequence );

    const auto read2 = std::make_shared< Read >( BasePairSequence( 4, 'A' ), std::string( 4, 'Q' ), "0", Cigar( "4M" ),
                                                 0, startPos + 1, 0, 0, 0, 0, 0, refSequence );

    echidna::io::readIntervalTree_t readContainer1( 0, 100 );
    readContainer1.insert( read1 );
    RegionsReads regionSetReads1( region, readContainer1.getFullRange(), 0 );

    echidna::io::readIntervalTree_t readContainer2( 0, 100 );
    readContainer2.insert( read2 );
    RegionsReads regionSetReads2( region, readContainer2.getFullRange(), 0 );

    echidna::io::perSampleRegionsReads_t perSampleReads = {{"sample1", regionSetReads1}, {"sample2", regionSetReads2}};

    const auto refCalls = buildRefCall( region, perSampleReads, 10, {2, 2}, 0.2 );

    BOOST_REQUIRE_EQUAL( refCalls.size(), 3 );

    BOOST_CHECK_EQUAL( refCalls[0].getAnnotation( Annotation::BEG ), 1 );
    BOOST_CHECK_EQUAL( refCalls[0].getAnnotation( Annotation::END ), 1 );
    BOOST_CHECK_EQUAL( refCalls[0].getAnnotation( Annotation::LEN ), 1 );
    BOOST_CHECK_EQUAL( refCalls[0].samples[0].getAnnotation( Annotation::MIN_DP ), 1 );
    BOOST_CHECK_EQUAL( refCalls[0].samples[0].getAnnotation( Annotation::FORMAT_DP ), 1 );
    BOOST_CHECK_CLOSE( refCalls[0].samples[0].getAnnotation( Annotation::GQ ), 3.0103, 1 );
    BOOST_CHECK_EQUAL( refCalls[0].samples[1].getAnnotation( Annotation::MIN_DP ), 0 );
    BOOST_CHECK_EQUAL( refCalls[0].samples[1].getAnnotation( Annotation::FORMAT_DP ), 0 );
    BOOST_CHECK_EQUAL( std::isnan( refCalls[0].samples[1].getAnnotation( Annotation::GQ ) ), true );

    BOOST_CHECK_EQUAL( refCalls[1].getAnnotation( Annotation::BEG ), 2 );
    BOOST_CHECK_EQUAL( refCalls[1].getAnnotation( Annotation::END ), 3 );
    BOOST_CHECK_EQUAL( refCalls[1].getAnnotation( Annotation::LEN ), 2 );
    BOOST_CHECK_EQUAL( refCalls[1].samples[0].getAnnotation( Annotation::MIN_DP ), 1 );
    BOOST_CHECK_EQUAL( refCalls[1].samples[0].getAnnotation( Annotation::FORMAT_DP ), 1 );
    BOOST_CHECK_CLOSE( refCalls[1].samples[0].getAnnotation( Annotation::GQ ), 3.0103, 1 );
    BOOST_CHECK_EQUAL( refCalls[1].samples[1].getAnnotation( Annotation::MIN_DP ), 1 );
    BOOST_CHECK_EQUAL( refCalls[1].samples[1].getAnnotation( Annotation::FORMAT_DP ), 1 );
    BOOST_CHECK_CLOSE( refCalls[1].samples[1].getAnnotation( Annotation::GQ ), 3.0103, 1 );

    BOOST_CHECK_EQUAL( refCalls[2].getAnnotation( Annotation::BEG ), 4 );
    BOOST_CHECK_EQUAL( refCalls[2].getAnnotation( Annotation::END ), 5 );
    BOOST_CHECK_EQUAL( refCalls[2].getAnnotation( Annotation::LEN ), 2 );
    BOOST_CHECK_EQUAL( refCalls[2].samples[0].getAnnotation( Annotation::MIN_DP ), 0 );
    BOOST_CHECK_EQUAL( refCalls[2].samples[0].getAnnotation( Annotation::FORMAT_DP ), 0 );
    BOOST_CHECK_EQUAL( std::isnan( refCalls[2].samples[0].getAnnotation( Annotation::GQ ) ), true );
    BOOST_CHECK_EQUAL( refCalls[2].samples[1].getAnnotation( Annotation::MIN_DP ), 1 );
    BOOST_CHECK_EQUAL( refCalls[2].samples[1].getAnnotation( Annotation::FORMAT_DP ), 1 );
    BOOST_CHECK_CLOSE( refCalls[2].samples[1].getAnnotation( Annotation::GQ ), 3.0103, 1 );
}
