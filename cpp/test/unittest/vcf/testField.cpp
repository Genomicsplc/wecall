// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <set>
#include "vcf/field.hpp"

BOOST_AUTO_TEST_CASE( testInfoFieldsContainNoDuplicates )
{
    const auto infoKeys = wecall::vcf::info::getVCFKeys( true );
    const std::set< std::string > asSet( infoKeys.cbegin(), infoKeys.cend() );
    BOOST_CHECK_EQUAL( asSet.size(), infoKeys.size() );
}

BOOST_AUTO_TEST_CASE( testCanConstructAVCFInfoFieldForEachKey )
{
    const auto infoKeys = wecall::vcf::info::getVCFKeys( true );
    for ( const auto & infoKey : infoKeys )
    {
        const auto formatedField = wecall::vcf::Field::infoFieldFromID( infoKey );
        BOOST_CHECK_EQUAL( formatedField.m_id, infoKey );
    }
}
