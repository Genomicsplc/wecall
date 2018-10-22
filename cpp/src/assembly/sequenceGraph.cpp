// All content Copyright (C) 2018 Genomics plc
#include "assembly/sequenceGraph.hpp"
#include "io/read.hpp"

#include <stdexcept>
#include <cmath>
#include <cassert>
#include <iostream>
#include <utility>
#include <queue>
#include <set>
#include <stack>
#include <boost/algorithm/string/join.hpp>

namespace wecall
{
namespace assembly
{

    std::size_t SequenceGraph::Chain::getSupport() const
    {
        if ( this->size() == 1 )
        {
            return 0;
        }
        else
        {
            std::size_t support = 0;
            for ( std::size_t index = 1; index < m_nodes.size(); ++index )
            {
                const char outEdge = m_nodes[index]->lastCharacter();
                const std::size_t edgeSupport = m_nodes[index - 1]->outSupport( outEdge );
                support = std::max( support, edgeSupport );
            }
            return support;
        }
    }

    bool SequenceGraph::KmerSequence::operator<( const KmerSequence & other ) const
    {
        if ( m_length != other.m_length )
        {
            return m_length < other.m_length;
        }
        for ( std::size_t i = 0; i < m_length; ++i )
        {
            if ( m_orginalSequence.at( m_pos + i ) != other.m_orginalSequence.at( other.m_pos + i ) )
            {
                return m_orginalSequence.at( m_pos + i ) < other.m_orginalSequence.at( other.m_pos + i );
            }
        }

        return false;
    }

    bool SequenceGraph::KmerSequence::operator==( const KmerSequence & other ) const
    {
        return m_length == other.m_length and
               std::equal( m_orginalSequence.cbegin() + m_pos, m_orginalSequence.cbegin() + m_pos + m_length,
                           other.m_orginalSequence.cbegin() + other.m_pos );
    }

    void SequenceGraph::Chain::computeSequence()
    {
        const auto kmerSize = m_initialSequence.size();
        auto sequence_str = m_initialSequence.str() + std::string( m_nodes.size() - 1, 'N' );

        for ( std::size_t index = 1; index < m_nodes.size(); ++index )
        {
            sequence_str[index - 1 + kmerSize] = m_nodes[index]->lastCharacter();
        }
        m_sequence = utils::BasePairSequence( sequence_str );
        m_endSequence = m_sequence.substr( m_sequence.size() - m_initialSequence.size(), m_initialSequence.size() );
    }

    SequenceGraph::Chain & SequenceGraph::Chain::operator+=( const Chain & other )
    {
        assert( not m_nodes.empty() );
        m_nodes.insert( m_nodes.end(), other.m_nodes.cbegin() + 1, other.m_nodes.cend() );
        m_sequence = m_sequence +
                     other.m_sequence.substr( other.m_initialSequence.size(),
                                              other.m_sequence.size() - other.m_initialSequence.size() );
        m_endSequence = other.m_endSequence;
        return *this;
    }

    std::vector< variant::varPtr_t > SequenceGraph::Chain::getVariants(
        const utils::referenceSequencePtr_t & refSeq ) const
    {
        std::vector< variant::varPtr_t > variantSet;
        const auto kmerSize = m_initialSequence.size();
        const auto startRefPositions = m_nodes.front()->getRefPositions();
        const auto endRefPositions = m_nodes.back()->getRefPositions();

        const auto alt = this->getSequence();

        for ( const std::size_t refStart : startRefPositions )
        {
            for ( const std::size_t refEnd : endRefPositions )
            {
                const std::size_t adjustedRefEnd = refEnd + kmerSize;
                if ( refStart <= adjustedRefEnd )
                {
                    const auto trimmed =
                        variant::Variant( refSeq, caller::Region( refSeq->contig(), refStart, adjustedRefEnd ), alt )
                            .getTrimmed();
                    if ( not trimmed->empty() )
                    {
                        const auto newVar = trimmed->getLeftAligned( refSeq->start() );
                        variantSet.push_back( newVar );
                    }
                    else
                    {
                        return {};
                    }
                }
            }
        }
        return variantSet;
    }

    SequenceGraph::SequenceGraph( const std::size_t kmerSize, const std::size_t minEdgeBaseQuality )
        : m_kmerSize( kmerSize ), m_minEdgeBaseQuality( minEdgeBaseQuality )
    {
    }

    //-----------------------------------------------------------------------------------------

    std::pair< std::shared_ptr< Node >, std::shared_ptr< Node > > SequenceGraph::addEdge(
        const utils::BasePairSequence * originalSequence,
        const std::size_t firstIndex,
        const std::size_t support )
    {
        /// Check if nodes already exist. If not, create them. Join node1 to node2 by an out-going
        /// edge from 1 to 2.

        const KmerSequence startNodeKmer( originalSequence, firstIndex, m_kmerSize );
        const KmerSequence endNodeKmer( originalSequence, firstIndex + 1, m_kmerSize );

        auto startNodeIt = m_nodeMap.find( startNodeKmer );
        auto endNodeIt = m_nodeMap.find( endNodeKmer );

        if ( startNodeIt == m_nodeMap.end() )
        {
            startNodeIt = m_nodeMap.emplace( startNodeKmer, std::make_shared< Node >( startNodeKmer.front(),
                                                                                      startNodeKmer.back() ) ).first;
        }

        if ( endNodeIt == m_nodeMap.end() )
        {
            endNodeIt = m_nodeMap.emplace( endNodeKmer,
                                           std::make_shared< Node >( endNodeKmer.front(), endNodeKmer.back() ) ).first;
        }

        const std::size_t adjustedSupport = ( support < m_minEdgeBaseQuality ) ? 0 : support;
        ( startNodeIt->second )->addOutEdge( endNodeIt->second, adjustedSupport );
        ( endNodeIt->second )->addInEdge( startNodeIt->second );
        return std::make_pair( startNodeIt->second, endNodeIt->second );
    }

    std::vector< SequenceGraph::Chain > SequenceGraph::getChains( const std::size_t minSupport ) const
    {
        std::vector< SequenceGraph::Chain > chains;
        const auto branchNodes = this->getBranchNodes();

        for ( const auto branchNode : branchNodes )
        {
            const auto firstSuccessors = branchNode.second->getSuccessors();
            for ( const auto & firstSuccessor : firstSuccessors )
            {
                // Always start with branch
                SequenceGraph::Chain newChain( branchNode.first.getSequence(), branchNode.second );

                auto currentNext = firstSuccessor;
                while ( currentNext->isRegular() )
                {
                    newChain.push_back( currentNext );
                    currentNext = currentNext->getSuccessors().front();
                }

                // Always end with branch.
                newChain.push_back( currentNext );

                if ( newChain.getSupport() >= minSupport )
                {
                    newChain.computeSequence();
                    chains.push_back( newChain );
                }
            }
        }

        return chains;
    }

    bool SequenceGraph::Chain::isAltSequence() const
    {
        if ( m_nodes.size() <= 2 and m_nodes.front()->isRef() and m_nodes.back()->isRef() )
        {
            return false;
        }

        if ( m_nodes.front()->isTerminal() or m_nodes.back()->isTerminal() )
        {
            return false;
        }

        for ( std::size_t index = 1; index < m_nodes.size() - 1; ++index )
        {
            if ( m_nodes[index]->isRef() )
            {
                return false;
            }
        }
        return true;
    }

    std::vector< SequenceGraph::Chain > SequenceGraph::getPathsBetweenRefNodes( const std::size_t minSupport ) const
    {
        typedef std::vector< std::size_t > path_t;
        const auto pureChains = this->getChains( minSupport );
        const auto hasSupport = [pureChains, minSupport]( const path_t & path ) -> bool
        {
            std::map< std::size_t, std::size_t > nRawChains;
            for ( const auto item : path )
            {
                ++nRawChains[item];
            }
            std::size_t support = 0;
            for ( const auto & pair : nRawChains )
            {
                support = std::max( support, pureChains[pair.first].getSupport() / pair.second );
            }
            return support >= minSupport;
        };

        std::map< utils::BasePairSequence, std::vector< std::size_t > > pureChainsFromNode;
        std::vector< utils::BasePairSequence > startNodes;

        std::vector< SequenceGraph::Chain > chainsForVariants;
        for ( std::size_t chainIndex = 0; chainIndex < pureChains.size(); ++chainIndex )
        {
            const auto & chain = pureChains[chainIndex];
            if ( chain.isAltSequence() )
            {
                pureChainsFromNode[chain.initialSequence()].push_back( chainIndex );
                if ( chain.front()->isRef() )
                {
                    startNodes.push_back( chain.initialSequence() );
                }
            }
            else if ( chain.size() >= 2 and chain.front()->isRef() and chain.back()->isRef() )
            {
                chainsForVariants.push_back( chain );
            }
        }

        for ( const auto startNode : startNodes )
        {
            std::map< utils::BasePairSequence, std::vector< path_t > > pathsToNode;

            std::priority_queue< utils::BasePairSequence > aQueueOfNodes;
            aQueueOfNodes.push( startNode );
            pathsToNode[startNode] = {{}};

            while ( not aQueueOfNodes.empty() )
            {
                const auto currentNode = aQueueOfNodes.top();
                aQueueOfNodes.pop();

                const auto pathsToCurrent = pathsToNode[currentNode];

                for ( const auto & chainsFromCurrent : pureChainsFromNode[currentNode] )
                {
                    const auto endNodeKmer = pureChains[chainsFromCurrent].endSequence();
                    auto & pathsToThisNode = pathsToNode[endNodeKmer];

                    if ( pathsToThisNode.size() <= 20 )
                    {
                        for ( const auto & pathToCurrent : pathsToCurrent )
                        {
                            auto newPath = pathToCurrent;
                            newPath.push_back( chainsFromCurrent );
                            if ( hasSupport( newPath ) )
                            {
                                pathsToThisNode.push_back( newPath );
                            }
                        }
                        aQueueOfNodes.push( endNodeKmer );
                    }
                }
            }

            for ( const auto & path : pathsToNode )
            {
                for ( const auto & chainIndicies : path.second )
                {
                    if ( not chainIndicies.empty() )
                    {
                        Chain chain = pureChains[chainIndicies[0]];
                        for ( std::size_t index = 1; index < chainIndicies.size(); ++index )
                        {
                            chain += pureChains[chainIndicies[index]];
                        }

                        if ( chain.size() >= 2 and chain.back()->isRef() )
                        {
                            chainsForVariants.push_back( chain );
                        }
                    }
                }
            }
        }
        return chainsForVariants;
    }

    SequenceGraph::nodeMap_t SequenceGraph::getBranchNodes() const
    {
        nodeMap_t branchNodes;
        for ( const auto & nodePair : m_nodeMap )
        {
            if ( nodePair.second->isBranch() )
            {
                branchNodes[nodePair.first] = nodePair.second;
            }
        }
        return branchNodes;
    }

    bool SequenceGraph::addReadSequence( const utils::BasePairSequence * readSequence,
                                         const utils::QualitySequence & qualitySequence,
                                         const bool disallowRepeats )
    {
        ECHIDNA_ASSERT( readSequence, "KmerSequence created from temporary object" );
        const std::size_t length = readSequence->size();
        assert( qualitySequence.size() == length );

        std::set< Node * > currentNodes;

        ECHIDNA_ASSERT( length > m_kmerSize, "Require sequence length to be longer than kmerSize" );
        bool hasRepeat = false;

        for ( std::size_t index = 0; index < length - m_kmerSize; ++index )
        {
            const auto edgeQual = static_cast< std::size_t >( qualitySequence[index + m_kmerSize] );
            const auto nodes = this->addEdge( readSequence, index, edgeQual );

            if ( currentNodes.count( nodes.first.get() ) > 0 or currentNodes.count( nodes.second.get() ) > 0 )
            {
                hasRepeat = true;
                if ( disallowRepeats )
                {
                    return hasRepeat;
                }
            }
            else
            {
                currentNodes.insert( nodes.first.get() );
            }
        }
        return hasRepeat;
    }

    bool SequenceGraph::labelAsReference( const std::shared_ptr< Node > & node, const int64_t refPos ) const
    {
        auto beforeReferencePositions = node->getRefPositions();
        node->addRefPos( refPos );
        if ( beforeReferencePositions.empty() )
        {
            return false;
        }
        else if ( beforeReferencePositions.count( refPos ) > 0 )
        {
            return beforeReferencePositions.size() != 1;
        }
        else
        {
            return true;
        }
    }

    //-----------------------------------------------------------------------------------------

    bool SequenceGraph::addReferenceSequence( const utils::referenceSequencePtr_t & referenceSequence,
                                              const caller::Region & region,
                                              const bool disallowRepeats )
    {
        const std::size_t length = region.size();

        bool repeatNode = false;

        ECHIDNA_ASSERT( referenceSequence->region().contains( region ), "" );
        ECHIDNA_ASSERT( length > m_kmerSize, "Require sequence length to be longer than kmerSize" );

        const std::size_t indentFromStart = int64_to_sizet( region.start() - referenceSequence->start() );

        for ( std::size_t i = indentFromStart; i < indentFromStart + length - m_kmerSize; ++i )
        {
            const auto nodes = this->addEdge( &referenceSequence->sequence(), i, 0 );

            if ( this->labelAsReference( nodes.first, referenceSequence->start() + i ) or
                 this->labelAsReference( nodes.second, referenceSequence->start() + i + 1 ) )
            {
                repeatNode = true;
                if ( disallowRepeats )
                {
                    return repeatNode;
                }
            }
        }
        return repeatNode;
    }

    //-----------------------------------------------------------------------------------------

    bool SequenceGraph::containsNode( const utils::BasePairSequence & theNode ) const
    {
        assert( theNode.size() == this->m_kmerSize );
        for ( const auto & pair : m_nodeMap )
        {
            if ( theNode == pair.first.getSequence() )
            {
                return true;
            }
        }
        return false;
    }

    std::string SequenceGraph::Chain::toString() const
    {
        std::stringstream sstr;
        sstr << m_initialSequence;
        if ( front()->isRef() )
        {
            sstr << " (Ref) ";
            for ( const auto refPos : front()->getRefPositions() )
            {
                sstr << " " << refPos;
            }
        }
        else
        {
            sstr << "       ";
        }
        sstr << " --> ";
        sstr << m_endSequence;
        if ( back()->isRef() )
        {
            sstr << " (Ref) ";
            for ( const auto refPos : back()->getRefPositions() )
            {
                sstr << " " << refPos;
            }
        }
        else
        {
            sstr << "       ";
        }
        sstr << " length: " << size();
        sstr << " support: " << getSupport();
        sstr << " isAltSequence " << isAltSequence();
        return sstr.str();
    }

    std::string SequenceGraph::toString() const
    {
        std::stringstream sstr;
        const auto chains = this->getChains( 1 );
        for ( const auto & chain : chains )
        {
            sstr << chain.toString() << std::endl;
        }
        return sstr.str();
    }

    //-----------------------------------------------------------------------------------------

    std::ostream & operator<<( std::ostream & stream, const Node & theNode )
    {
        stream << theNode.firstCharacter() << " ---> " << theNode.lastCharacter();
        return stream;
    }

    //-----------------------------------------------------------------------------------------
}
}
