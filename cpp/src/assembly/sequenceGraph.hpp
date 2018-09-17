// All content Copyright (C) 2018 Genomics plc
#ifndef SEQUENCE_GRAPH_HPP
#define SEQUENCE_GRAPH_HPP

#include <string>
#include <memory>
#include <vector>
#include <array>
#include <map>
#include <utils/referenceSequence.hpp>
#include "assembly/node.hpp"
#include "utils/sequence.hpp"
#include "variant/type/variant.hpp"

namespace echidna
{

namespace assembly
{
    class SequenceGraph
    {
    public:
        class KmerSequence
            : boost::partially_ordered< KmerSequence, KmerSequence, boost::equality_comparable< KmerSequence > >
        {
        public:
            KmerSequence( const utils::BasePairSequence * origionalSequence,
                          const std::size_t pos,
                          const std::size_t length )
                : m_orginalSequence( *origionalSequence ), m_pos( pos ), m_length( length )
            {
                ECHIDNA_ASSERT( origionalSequence, "KmerSequence created from temporary object" );
            }

            bool operator<( const KmerSequence & other ) const;
            bool operator==( const KmerSequence & other ) const;

            utils::BasePairSequence getSequence() const { return m_orginalSequence.substr( m_pos, m_length ); }
            char front() const { return m_orginalSequence.at( m_pos ); }
            char back() const { return m_orginalSequence.at( m_pos + m_length - 1 ); }

        private:
            const utils::BasePairSequence & m_orginalSequence;
            const std::size_t m_pos;
            const std::size_t m_length;
        };

        typedef std::map< KmerSequence, std::shared_ptr< Node > > nodeMap_t;

        class Chain
        {
        public:
            Chain( const utils::BasePairSequence & initialSequence, const std::shared_ptr< Node > & initialNode )
                : m_initialSequence( initialSequence ), m_nodes( {initialNode} )
            {
            }
            Chain & operator+=( const Chain & other );

            std::string toString() const;

            std::vector< variant::varPtr_t > getVariants( const utils::referenceSequencePtr_t & refSeq ) const;
            const utils::BasePairSequence & initialSequence() const { return m_initialSequence; }
            const utils::BasePairSequence & endSequence() const { return m_endSequence; }

            const utils::BasePairSequence & getSequence() const { return m_sequence; }
            void computeSequence();

            void push_back( const std::shared_ptr< Node > & nextNode ) { m_nodes.push_back( nextNode ); }
            std::size_t getSupport() const;
            bool isAltSequence() const;

            const std::shared_ptr< Node > & front() const { return m_nodes.front(); }
            const std::shared_ptr< Node > & back() const { return m_nodes.back(); }
            std::size_t size() const { return m_nodes.size(); }

        private:
            const utils::BasePairSequence m_initialSequence;
            utils::BasePairSequence m_sequence;
            utils::BasePairSequence m_endSequence;
            std::vector< std::shared_ptr< Node > > m_nodes;
        };

    public:
        /// Constructor
        ///
        /// @kmerSize The size of kmer to use for Nodes in this graph. TODO: Should be able to deal with multiple kmer
        /// sizes?
        SequenceGraph( const std::size_t kmerSize, const std::size_t minEdgeBaseQuality );

        /// Destructor
        ~SequenceGraph() {}

        // Chain: A linearly linked walk of nodes between two branch nodes.
        // Chains includes branch nodes at each end but no other branch nodes.
        // Use to find bubbles. A bubble arises from a chain with non-terminal branch nodes at either end.
        // Use terminal nodes to find strands of sequence arising from large cnv breakpoints.
        std::vector< SequenceGraph::Chain > getChains( const std::size_t minSupport ) const;

        std::vector< SequenceGraph::Chain > getPathsBetweenRefNodes( const std::size_t minSupport ) const;

        // Warning: Requires readSequence to not go out of scope for duration of sequence graph.
        bool addReadSequence( const utils::BasePairSequence * readSequence,
                              const utils::QualitySequence & qualitySequence,
                              const bool disallowRepeats );

        /// Add a reference sequence to the graph. All kmers will be added.
        // So, a node would contain a set of reference positions.
        bool addReferenceSequence( const utils::referenceSequencePtr_t & referenceSequence,
                                   const caller::Region & region,
                                   const bool disallowRepeats );

        /// Check if this graph already contains a specified Node.
        /// @param theNode Node to check
        /// @return true if present otherwise false
        bool containsNode( const utils::BasePairSequence & theNode ) const;

        std::string toString() const;

        std::size_t size() const { return m_nodeMap.size(); }
        std::size_t kmerSize() const { return m_kmerSize; }

    private:
        nodeMap_t getBranchNodes() const;

        bool labelAsReference( const std::shared_ptr< Node > & node, const int64_t refPos ) const;

        std::pair< std::shared_ptr< Node >, std::shared_ptr< Node > > addEdge(
            const utils::BasePairSequence * originalSequence,
            const std::size_t firstIndex,
            const std::size_t support );

        const std::size_t m_kmerSize;
        const std::size_t m_minEdgeBaseQuality;

        nodeMap_t m_nodeMap;
    };

    /// Overloaded iostream operators
    std::ostream & operator<<( std::ostream & stream, const Node & theNode );
}
}

#endif
