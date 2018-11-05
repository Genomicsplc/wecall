// All content Copyright (C) 2018 Genomics plc
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "assembly/node.hpp"
#include "common.hpp"

using wecall::assembly::Node;

BOOST_AUTO_TEST_CASE( testIsolatedNode )
{
    auto node = std::make_shared< Node >( 'A', 'C' );

    BOOST_CHECK_EQUAL( node->lastCharacter(), 'C' );
    BOOST_CHECK( not node->isRegular() );
    BOOST_CHECK( node->isBranch() );
    BOOST_CHECK( node->isTerminal() );
}

BOOST_AUTO_TEST_CASE( testRegularNode )
{
    auto node = std::make_shared< Node >( 'A', 'C' );
    auto inNode = std::make_shared< Node >( 'G', 'T' );
    auto outNode = std::make_shared< Node >( 'T', 'G' );

    node->addInEdge( inNode );
    node->addOutEdge( outNode, 0 );

    BOOST_CHECK( node->isRegular() );
    BOOST_CHECK( not node->isBranch() );
    BOOST_CHECK( not node->isTerminal() );
    BOOST_CHECK_EQUAL( ( *node->begin() ).lastCharacter(), 'G' );
}
