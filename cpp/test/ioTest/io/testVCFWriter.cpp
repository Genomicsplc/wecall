#include "io/vcfWriter.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <deque>
#include <vcf/reader.hpp>
#include "variant/type/variant.hpp"
#include "utils/date.hpp"
#include "ioFixture.hpp"
#include "unittest/vcf/VCFTestUtils.hpp"

std::string readEntireFileIntoStdString( std::string filename )
{
    std::ifstream fp( filename );
    std::stringstream sstr;
    sstr << fp.rdbuf();
    return sstr.str();
}

// BOOST_FIXTURE_TEST_CASE( testCantWriteHeaderTwice, echidna::test::FastaFileFixture )
//{
//    std::string tempFileTemplate = "/tmp/tmp_vcf_output_XXXXXX";
//    const std::string tempFilename = tempFileTemplate.c_str();
//
//    echidna::io::VCFWriter writer( tempFilename, false, false, false );
//
//    echidna::caller::params::Application applicationParams( "Edna", "1.0", "blah", "25-11-1987", "options" );
//    std::vector< std::string > sampleNames{"NA17287", "NA17291"};
//    std::vector< echidna::vcf::FilterDesc > filterDescs{{"SB", "this is a test filter"}};
//    std::vector< echidna::caller::Region > mapContigsToLengths;
//    mapContigsToLengths.emplace_back( "2", 1, 20 );
//    mapContigsToLengths.emplace_back( "10", 1, 100 );
//    mapContigsToLengths.emplace_back( "20", 1, 20 );
//    writer.writeHeader( "VCF4.2", applicationParams, refFilename, sampleNames, filterDescs, mapContigsToLengths );
//
//    BOOST_CHECK_THROW(
//        writer.writeHeader( "VCF4.2", applicationParams, refFilename, sampleNames, filterDescs, mapContigsToLengths ),
//        echidna::utils::echidna_exception );
//}
//
// BOOST_AUTO_TEST_CASE( compileInfoShouldCombineConcatenationAnnotations )
//{
//    echidna::variant::varPtr_t variant1( std::make_shared< echidna::variant::SNP >( "20", 66370, "G", "A" ) );
//    echidna::variant::varPtr_t variant2( std::make_shared< echidna::variant::SNP >( "20", 66370, "G", "T" ) );
//
//    echidna::caller::CallSet callset( 1 );
//    auto & call1 = callset.addVarCall( variant1, 101.0, {{echidna::caller::Call::Type::VAR}} );
//    call1.addAnnotation( echidna::Annotation::VC, 22l );
//    auto & call2 = callset.addVarCall( variant2, 102.0, {{echidna::caller::Call::Type::VAR}} );
//    call2.addAnnotation( echidna::Annotation::VC, 44l );
//
//    echidna::vcf::Info procesedInfo = echidna::io::VCFWriter::compileInfo( callset.begin(), callset.end() );
//    echidna::vcf::Info expectedInfo = {
//        {"VC", {"22", "44"}},
//    };
//    check_equal_info( expectedInfo, procesedInfo );
//}
//
// BOOST_AUTO_TEST_CASE( compileInfoShouldCombineFirstAnnotations )
//{
//    echidna::variant::varPtr_t variant1( std::make_shared< echidna::variant::SNP >( "20", 66370, "G", "A" ) );
//    echidna::variant::varPtr_t variant2( std::make_shared< echidna::variant::SNP >( "20", 66370, "G", "T" ) );
//
//    echidna::caller::CallSet callset( 1 );
//    auto & call1 = callset.addVarCall( variant1, 101.0, {{echidna::caller::Call::Type::VAR}} );
//    call1.addAnnotation( echidna::Annotation::DP, 1l );
//    auto & call2 = callset.addVarCall( variant2, 102.0, {{echidna::caller::Call::Type::VAR}} );
//    call2.addAnnotation( echidna::Annotation::DP, 2l );
//
//    echidna::vcf::Info procesedInfo = echidna::io::VCFWriter::compileInfo( callset.begin(), callset.end() );
//    echidna::vcf::Info expectedInfo = {
//        {"DP", {"1"}},
//    };
//    check_equal_info( expectedInfo, procesedInfo );
//}
//
// BOOST_FIXTURE_TEST_CASE( idFieldShouldBeWrittenCorrectlyForSNPs, echidna::test::FastaFileFixture )
//{
//    using namespace echidna::vcf;
//    std::string tempFileTemplate = "/tmp/tmp_vcf_output_XXXXXX";
//    const std::string tempFilename = tempFileTemplate.c_str();
//
//    echidna::io::VCFWriter writer( tempFilename, false, false, false );
//
//    echidna::caller::params::Application applicationParams( "Edna", "1.0", "blah", "25-11-1987", "options" );
//    std::vector< std::string > sampleNames{"NA17287", "NA17291"};
//    std::vector< echidna::vcf::FilterDesc > filterDescs{{"REFCALL", ""},
//                                                        {"FOO", "A filter for testing that filter data is written"}};
//    writer.writeHeader( "VCF4.2", applicationParams, refFilename, sampleNames, filterDescs, {} );
//
//    echidna::variant::varPtr_t snpWithoutVarIds( std::make_shared< echidna::variant::SNP >( "20", 66370, "G", "A" ) );
//    echidna::variant::varPtr_t snpWithVarIds( std::make_shared< echidna::variant::SNP >( "20", 72132, "A", "G" ) );
//
//    echidna::caller::CallSet callset( 1 );
//    callset.addRefCall( 66370, 42.0, {{echidna::caller::Call::Type::REF}} );
//    callset.addVarCall( snpWithoutVarIds, 101.0, {{echidna::caller::Call::Type::VAR}} );
//    std::set< std::string > testIds{"testId1", "testId2"};
//    callset.addVarCall( snpWithVarIds, 102.0, {{echidna::caller::Call::Type::VAR}} ).varIds = testIds;
//
//    writer.contig( "20" );
//    writer.writeCallSet( ref, callset );
//
//    echidna::vcf::Reader reader( tempFilename );
//
//    std::vector< std::set< std::string > > expectedVarIds{{}, {}, testIds};
//    std::vector< std::set< std::string > > obtainedVarIds;
//    while ( auto record = reader.getNextRecord() )
//    {
//        obtainedVarIds.push_back( record.get().m_ids );
//    }
//    BOOST_CHECK_EQUAL( expectedVarIds.size(), obtainedVarIds.size() );
//    for ( auto i = 0; i < expectedVarIds.size(); ++i )
//    {
//        BOOST_CHECK( expectedVarIds[i] == obtainedVarIds[i] );
//    }
//}
//
// BOOST_FIXTURE_TEST_CASE( idFieldShouldBeWrittenCorrectlyForINSs, echidna::test::FastaFileFixture )
//{
//    using namespace echidna::vcf;
//    std::string tempFileTemplate = "/tmp/tmp_vcf_output_XXXXXX";
//    const std::string tempFilename = tempFileTemplate.c_str();
//
//    echidna::io::VCFWriter writer( tempFilename, false, false, false );
//
//    echidna::caller::params::Application applicationParams( "Edna", "1.0", "blah", "25-11-1987", "options" );
//    std::vector< std::string > sampleNames{"NA17287", "NA17291"};
//    std::vector< echidna::vcf::FilterDesc > filterDescs{{"REFCALL", ""},
//                                                        {"FOO", "A filter for testing that filter data is written"}};
//    writer.writeHeader( "VCF4.2", applicationParams, refFilename, sampleNames, filterDescs, {} );
//
//    echidna::variant::varPtr_t snpWithoutVarIds( std::make_shared< echidna::variant::SNP >( "20", 66370, "G", "A" ) );
//    echidna::variant::varPtr_t insWithVarIds( std::make_shared< echidna::variant::Insertion >( "20", 72133, "TTG" ) );
//
//    echidna::caller::CallSet callset( 1 );
//    callset.addRefCall( 66370, 42.0, {{echidna::caller::Call::Type::REF}} );
//    callset.addVarCall( snpWithoutVarIds, 101.0, {{echidna::caller::Call::Type::VAR}} );
//    std::set< std::string > testIds{"testId1", "testId2"};
//    callset.addVarCall( insWithVarIds, 102.0, {{echidna::caller::Call::Type::VAR}} ).varIds = testIds;
//
//    writer.contig( "20" );
//    writer.writeCallSet( ref, callset );
//
//    echidna::vcf::Reader reader( tempFilename );
//
//    std::vector< std::set< std::string > > expectedVarIds{{}, {}, testIds};
//    std::vector< std::set< std::string > > obtainedVarIds;
//    while ( auto record = reader.getNextRecord() )
//    {
//        obtainedVarIds.push_back( record.get().m_ids );
//    }
//    BOOST_CHECK_EQUAL( expectedVarIds.size(), obtainedVarIds.size() );
//    for ( auto i = 0; i < expectedVarIds.size(); ++i )
//    {
//        BOOST_CHECK( expectedVarIds[i] == obtainedVarIds[i] );
//    }
//}
