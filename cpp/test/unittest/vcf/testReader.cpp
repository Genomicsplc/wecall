#define BOOST_TEST_DYN_LINK
#include <cstdlib>
#include <tuple>
#include <boost/test/unit_test.hpp>
#include <boost/algorithm/string.hpp>

#include "vcf/reader.hpp"
#include "io/tabixVCFFile.hpp"
#include "vcf/filterDescription.hpp"
#include "VCFTestUtils.hpp"

using std::string;

BOOST_AUTO_TEST_CASE( shouldParseValidVCFFilterHeaderUpperCaseID )
{
    const std::string line = "##FILTER=<ID=ABCDEFGHIJKLMNOPQRSTUVWXYZ,Description=\"A description of a filter.\">";
    echidna::vcf::FilterDesc filterDesc = echidna::io::TabixVCFFile::parseFilterHeaderLine( line );
    BOOST_CHECK_EQUAL( filterDesc.id, "ABCDEFGHIJKLMNOPQRSTUVWXYZ" );
    BOOST_CHECK_EQUAL( filterDesc.description, "A description of a filter." );
}

BOOST_AUTO_TEST_CASE( shouldParseValidVCFFilterHeaderLowerCaseID )
{
    const std::string line = "##FILTER=<ID=abcdefghijklmnopqrstuvwxyz,Description=\"A description of a filter.\">";
    echidna::vcf::FilterDesc filterDesc = echidna::io::TabixVCFFile::parseFilterHeaderLine( line );
    BOOST_CHECK_EQUAL( filterDesc.id, "abcdefghijklmnopqrstuvwxyz" );
    BOOST_CHECK_EQUAL( filterDesc.description, "A description of a filter." );
}

BOOST_AUTO_TEST_CASE( shouldParseValidVCFFilterHeaderDigitsAndPunctuation )
{
    const std::string line = "##FILTER=<ID=_1234567890,Description=\"A description of a filter.\">";
    echidna::vcf::FilterDesc filterDesc = echidna::io::TabixVCFFile::parseFilterHeaderLine( line );
    BOOST_CHECK_EQUAL( filterDesc.id, "_1234567890" );
    BOOST_CHECK_EQUAL( filterDesc.description, "A description of a filter." );
}

BOOST_AUTO_TEST_CASE( shouldRaiseOnINFOHeaderType )
{
    const std::string line = "##INFO=<ID=AA,Number=A,Type=String,Description=\"Genotype\">";
    BOOST_CHECK_THROW( echidna::vcf::FilterDesc filterDesc = echidna::io::TabixVCFFile::parseFilterHeaderLine( line ),
                       echidna::utils::echidna_exception );
}

BOOST_AUTO_TEST_CASE( shouldRaiseOnFORMATHeaderType )
{
    const std::string line = "##FORMAT=<ID=AA,Number=A,Type=String,Description=\"Genotype\">";
    BOOST_CHECK_THROW( echidna::vcf::FilterDesc filterDesc = echidna::io::TabixVCFFile::parseFilterHeaderLine( line ),
                       echidna::utils::echidna_exception );
}

BOOST_AUTO_TEST_CASE( shouldRaiseOnInvalidVCFFilterHeader )
{
    const std::string line = "##FILTER=<ID=q10;Description=\"A description of a filter.\">";
    BOOST_CHECK_THROW( echidna::vcf::FilterDesc filterDesc = echidna::io::TabixVCFFile::parseFilterHeaderLine( line ),
                       echidna::utils::echidna_exception );
}

BOOST_AUTO_TEST_CASE( shouldDetectFilterIdInFilterDescs )
{
    using namespace echidna::vcf;

    std::string id1 = "gq1";
    std::string id2 = "gq2";
    std::set< FilterDesc > filterDescs = {FilterDesc( id1, "Random description" )};

    BOOST_CHECK( echidna::io::TabixVCFFile::containsFilterId( filterDescs, id1 ) );
    BOOST_CHECK( not echidna::io::TabixVCFFile::containsFilterId( filterDescs, id2 ) );
}
