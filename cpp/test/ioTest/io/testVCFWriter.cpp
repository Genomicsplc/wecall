// All content Copyright (C) 2018 Genomics plc
#include "io/vcfWriter.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <deque>
#include <vcf/reader.hpp>
#include <io/tabixVCFFile.hpp>
#include "variant/type/variant.hpp"
#include "utils/date.hpp"
#include "ioFixture.hpp"

using echidna::variant::Variant;
using echidna::utils::ReferenceSequence;
using echidna::caller::Region;
using echidna::caller::Call;
using echidna::caller::Annotation;

std::string readEntireFileIntoStdString( const std::string & filename )
{
    std::ifstream fp( filename );
    std::stringstream sstr;
    sstr << fp.rdbuf();
    return sstr.str();
}

void check_equal_info( const echidna::vcf::Info & expected, const echidna::vcf::Info & parsed )
{
    std::ostringstream info_size_message;
    info_size_message << "Different length info fields (expected " << expected.size() << " elements: [";
    for ( auto & item : expected )
    {
        info_size_message << item.first << ", ";
    }
    info_size_message << "]; got " << parsed.size() << " elements: [";
    for ( auto & item : parsed )
    {
        info_size_message << item.first << ", ";
    }
    info_size_message << "]).";
    BOOST_CHECK_MESSAGE( expected.size() == parsed.size(), info_size_message.str().c_str() );
    for ( unsigned int i = 0; i != std::min( expected.size(), parsed.size() ); ++i )
    {
        std::ostringstream item_key_message;
        item_key_message << "Different keys (expected '" << expected[i].first << "'; got '" << parsed[i].first
                         << "') for i = " << i;
        BOOST_CHECK_MESSAGE( expected[i].first == parsed[i].first, item_key_message.str().c_str() );
        std::ostringstream inner_length_message;
        inner_length_message << "Inner vectors of different size (expected '" << expected[i].second.size() << "'; got '"
                             << parsed[i].second.size() << "') for i = " << i;
        BOOST_CHECK_MESSAGE( expected[i].second.size() == parsed[i].second.size(), inner_length_message.str().c_str() );
        for ( unsigned int j = 0; j != std::min( expected[i].second.size(), parsed[i].second.size() ); ++j )
        {
            std::ostringstream inner_element_message;
            inner_element_message << "Elements don't match (expected '" << expected[i].second[j] << "'; got '"
                                  << parsed[i].second[j] << "') for i, j = " << i << ", " << j;
            BOOST_CHECK_MESSAGE( expected[i].second[j] == parsed[i].second[j], inner_element_message.str().c_str() );
        }
    }
}

BOOST_FIXTURE_TEST_CASE( testCantWriteHeaderTwice, echidna::test::FastaFileFixture )
{
    const std::string tempFilename = boost::filesystem::temp_directory_path().generic_string() + "/vcf_test";
    echidna::io::VCFWriter writer( tempFilename, false, false );

    echidna::caller::params::Application applicationParams( "Edna", "1.0", "blah", "25-11-1987", "options" );
    std::vector< std::string > sampleNames{"NA17287", "NA17291"};
    std::vector< echidna::vcf::FilterDesc > filterDescs{{"SB", "this is a test filter"}};
    std::vector< echidna::caller::Region > mapContigsToLengths;
    mapContigsToLengths.emplace_back( "2", 1, 20 );
    mapContigsToLengths.emplace_back( "10", 1, 100 );
    mapContigsToLengths.emplace_back( "20", 1, 20 );
    writer.writeHeader( "VCF4.2", applicationParams, refFilename, sampleNames, filterDescs, mapContigsToLengths );

    BOOST_CHECK_THROW(
        writer.writeHeader( "VCF4.2", applicationParams, refFilename, sampleNames, filterDescs, mapContigsToLengths ),
        echidna::utils::echidna_exception );
}

BOOST_AUTO_TEST_CASE( compileInfoForSingleVariant )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'G' ) );
    echidna::variant::varPtr_t variant = std::make_shared< Variant >( refSequence, Region( "1", 1, 2 ), "A" );

    Call call( variant, variant->interval(), 101.0, 1, {{Call::REF, Call::VAR}} );
    call.addAnnotation( Annotation::VC, 22l );
    echidna::caller::callVector_t calls = {call};

    echidna::vcf::Info processedInfo = echidna::io::VCFWriter::compileInfo( calls.cbegin() );
    echidna::vcf::Info expectedInfo = {{"VC", {"22"}}};
    check_equal_info( expectedInfo, processedInfo );
}

BOOST_AUTO_TEST_CASE( compileTwoInfoFieldsForSingleVariant )
{
    auto refSequence =
        std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 0, 5 ), std::string( 5, 'G' ) );
    echidna::variant::varPtr_t variant = std::make_shared< Variant >( refSequence, Region( "1", 1, 2 ), "A" );

    Call call( variant, variant->interval(), 101.0, 1, {{Call::REF, Call::VAR}} );
    call.addAnnotation( Annotation::VC, 22l );
    call.addAnnotation( Annotation::DP, 1l );
    echidna::caller::callVector_t calls = {call};

    echidna::vcf::Info procesedInfo = echidna::io::VCFWriter::compileInfo( calls.cbegin() );
    echidna::vcf::Info expectedInfo = {{"VC", {"22"}}, {"DP", {"1"}}};
    check_equal_info( expectedInfo, procesedInfo );
}

BOOST_FIXTURE_TEST_CASE( idFieldShouldBeWrittenCorrectlyForSNPs, echidna::test::FastaFileFixture )
{
    using namespace echidna::vcf;

    const std::string tempFilename = boost::filesystem::temp_directory_path().generic_string() + "/vcf_test";
    echidna::io::VCFWriter writer( tempFilename, false, false );

    echidna::caller::params::Application applicationParams( "Edna", "1.0", "blah", "25-11-1987", "options" );
    std::vector< std::string > sampleNames{"NA17287", "NA17291"};
    std::vector< echidna::vcf::FilterDesc > filterDescs{{"REFCALL", ""},
                                                        {"FOO", "A filter for testing that filter data is written"}};
    writer.writeHeader( "VCF4.2", applicationParams, refFilename, sampleNames, filterDescs, {} );

    auto refSequence = std::make_shared< echidna::utils::ReferenceSequence >( Region( "1", 60, 65 ), "CGCAG" );
    echidna::variant::varPtr_t variant = std::make_shared< Variant >( refSequence, Region( "1", 60, 61 ), "A" );

    Call call( variant, variant->interval(), 101.0, 1, {{Call::REF, Call::VAR}} );
    call.addAnnotation( Annotation::VC, 22l );
    echidna::caller::callVector_t calls = {call};

    writer.contig( "1" );
    refFiles.emplace_back( new echidna::io::FastaFile( refFilename ) );
    writer.writeCallSet( *( refFiles.back() ), calls );

    std::string s = readEntireFileIntoStdString( tempFilename );
    BOOST_CHECK( s.find( "1\t61\t.\tC\tA\t101\tPASS\tVC=22\tGT\t0/1" ) != s.size() );
}
