// All content Copyright (C) 2018 Genomics plc
#include <memory>
#include <string>
#include <map>

#include "assembly/node.hpp"
#include <cassert>

namespace echidna
{
namespace assembly
{

    void Node::addInEdge( const std::shared_ptr< Node > & inEdge )
    {
        const auto priorCharacter = inEdge->firstCharacter();
        m_inEdges.insert( priorCharacter );
    }

    //-----------------------------------------------------------------------------------------

    void Node::addOutEdge( const std::shared_ptr< Node > & outEdge, const std::size_t support )
    {
        const auto & nextCharacter = outEdge->lastCharacter();
        m_outEdges[nextCharacter] = outEdge;
        m_outCount[nextCharacter] += support;
    }

    std::vector< std::shared_ptr< Node > > Node::getSuccessors() const
    {
        std::vector< std::shared_ptr< Node > > successors;
        successors.reserve( m_outEdges.size() );
        for ( const auto & pair : m_outEdges )
        {
            successors.push_back( pair.second.lock() );
        }
        return successors;
    }
}
}
