// All content Copyright (C) 2018 Genomics plc
#ifndef SEQUENCE_NODE_HPP
#define SEQUENCE_NODE_HPP

#include <memory>
#include <string>
#include <map>
#include <set>
#include <vector>

#include "utils/sequence.hpp"

namespace echidna
{

namespace assembly
{
    /// A somwhat bloated class that represents a graph of k-mers (short DNA sequences). This is mainly meant
    /// for use in Echidna for assembling haplotypes from read data.
    const std::map< char, std::size_t > initCounts = {{'A', 0}, {'T', 0}, {'C', 0}, {'G', 0}};

    class Node
    {
    public:
        typedef std::map< char, std::weak_ptr< Node > > edges_t;

        class const_iterator
        {
        public:
            typedef Node value_type;

            const_iterator( edges_t::const_iterator it ) : m_it( it ) {}
            ~const_iterator() {}

            // TODO(ES): This is not safe as might have gone out of scope.
            const value_type & operator*() const { return *( m_it->second.lock() ); }
            const_iterator & operator++()
            {
                ++m_it;
                return *this;
            }
            bool operator!=( const_iterator rhs ) { return this->m_it != rhs.m_it; }
            bool operator==( const_iterator rhs ) { return this->m_it == rhs.m_it; }

        private:
            edges_t::const_iterator m_it;
        };

    public:
        /// Construct from BasePairSequence
        Node( const char firstCharacter, const char lastCharacter )
            : m_firstCharacter( firstCharacter ),
              m_lastCharacter( lastCharacter ),
              m_outCount( initCounts ),
              m_refPositions()
        {
        }

        /// Construct from Node
        Node( const Node & rhs )
            : m_firstCharacter( rhs.m_firstCharacter ),
              m_lastCharacter( rhs.m_lastCharacter ),
              m_inEdges( rhs.m_inEdges ),
              m_outEdges( rhs.m_outEdges ),
              m_outCount( rhs.m_outCount ),
              m_refPositions( rhs.m_refPositions )
        {
        }

        /// Move from Node
        Node( Node && rhs )
            : m_firstCharacter( std::move( rhs.m_firstCharacter ) ),
              m_lastCharacter( std::move( rhs.m_lastCharacter ) ),
              m_inEdges( std::move( rhs.m_inEdges ) ),
              m_outEdges( std::move( rhs.m_outEdges ) ),
              m_outCount( std::move( rhs.m_outCount ) ),
              m_refPositions( std::move( rhs.m_refPositions ) )
        {
        }

        /// Destructor
        ~Node() {}

        void addRefPos( const int64_t refPos ) { m_refPositions.insert( refPos ); };
        std::set< int64_t > getRefPositions() const { return m_refPositions; }
        bool isRef() const { return not m_refPositions.empty(); }

        void addInEdge( const std::shared_ptr< Node > & outEdge );
        void addOutEdge( const std::shared_ptr< Node > & outEdge, const std::size_t addSupport );

        std::size_t nPredecessors() const { return m_inEdges.size(); }
        std::size_t nSuccessors() const { return m_outEdges.size(); }

        std::size_t outSupport( const char outEdge ) { return m_outCount[outEdge]; }

        std::vector< std::shared_ptr< Node > > getSuccessors() const;

        bool isIsolated() const { return nPredecessors() == 0 and nSuccessors() == 0; }
        bool isRegular() const { return nPredecessors() == 1 and nSuccessors() == 1; }
        bool isBranch() const { return not isRegular(); }
        bool isTerminal() const { return m_inEdges.empty() or m_outEdges.empty(); }

        char firstCharacter() const { return m_firstCharacter; }
        char lastCharacter() const { return m_lastCharacter; }

        /// Return iterator to the first out edge
        const_iterator begin() const { return const_iterator( m_outEdges.begin() ); }

        /// Return iterator to the last+1 out edge
        const_iterator end() const { return const_iterator( m_outEdges.end() ); }

        /// Overloaded output operator
        friend std::ostream & operator<<( std::ostream & stream, const Node & theNode );

    private:
        const char m_firstCharacter;
        const char m_lastCharacter;

        std::set< char > m_inEdges;
        std::map< char, std::weak_ptr< Node > > m_outEdges;

        std::map< char, std::size_t > m_outCount;

        std::set< int64_t > m_refPositions;
    };
}
}
#endif
