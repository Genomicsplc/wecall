// All content Copyright (C) 2018 Genomics plc
#ifndef INTERVAL_TREE_HPP
#define INTERVAL_TREE_HPP

#include <functional>
#include <cassert>
#include <cstdint>
#include <memory>
#include <set>
#include <sstream>

#include "utils/logging.hpp"
#include "interval.hpp"
#include "utils/exceptions.hpp"
#include <boost/container/flat_set.hpp>

namespace wecall
{
namespace utils
{

    /// A binary tree that holds half-interval classes. Each node of the tree has a central point and stores (shared
    /// pointers to)
    // interval classes which contain that central point.
    /// Each node has a left and right bound. All intervals strickly contained in this range are either stored at this
    /// node
    // or within its children.
    template < typename _IntervalType,
               typename _IntervalLeftBoundaryComp,  ///<  These sort the intervals on each node.
               typename _IntervalRightBoundaryComp,
               typename _IntervalLeftBoundaryValue = _IntervalLeftBoundaryComp,  // Ideally should use default.
               typename _IntervalRightBoundaryValue = _IntervalRightBoundaryComp >
    class IntervalTree
    {
    public:
        class iterator;
        typedef std::pair< iterator, iterator > iterator_pair;
        typedef _IntervalType interval_type;
        typedef std::shared_ptr< interval_type > pointer;
        typedef interval_type & reference;

    private:
        class IntervalTreeNode;
        typedef std::shared_ptr< IntervalTreeNode > node_shared_ptr;
        typedef IntervalTreeNode * node_ptr;

        typedef _IntervalLeftBoundaryValue left_interval_value;
        typedef _IntervalRightBoundaryValue right_interval_value;
        typedef boost::container::flat_multiset< pointer, _IntervalLeftBoundaryComp > left_sorted_set;
        typedef boost::container::flat_multiset< pointer, _IntervalRightBoundaryComp > right_sorted_set;

        typedef typename boost::container::flat_multiset< pointer >::const_iterator
            set_iterator;  // Used for both left_sorted_set and right_sorted_set

    public:
        /// Construct an interval Tree for easily retrieving intervals that overlap a given interval
        ///
        /// @param left Default lower bound the interval tree spans. Gets decreased if necessary
        /// @param right Default upper bound the interval tree spans. Gets increased if necessary.
        IntervalTree( const int64_t & left, const int64_t & right )
            : m_root( std::make_shared< IntervalTreeNode >( Interval( left, right ), nullptr ) )
        {
        }

        /// Disabled copy constructor. No copying of interval trees.
        IntervalTree( const IntervalTree & rhs ) = delete;

        /// Disabled assigment operator. No copying of interval trees.
        IntervalTree & operator=( const IntervalTree & rhs ) = delete;

    public:
        /// Inserts a new interval pointer into the tree.
        /// interval will get inserted at the first node which has central value contained in the interval.
        ///
        /// @param intervalPtr Insert a std::shared_ptr to class of interval type.
        void insert( pointer intervalPtr )
        {
            const Interval newInterval( left_interval_value()( intervalPtr ), right_interval_value()( intervalPtr ) );
            this->widenTreeBoundary( newInterval );
            m_root->insert( intervalPtr, newInterval );
        }

        /// @return Returns true if IntervalTree is empty.
        bool empty() const { return m_root->empty(); }

        /// @return A string that identifies the tree for debugging purposes.
        std::string toString() const
        {
            std::stringstream sstrThisID;
            sstrThisID << "IntervalTree spanning region:- " << m_root->getBoundaries().toString();
            return sstrThisID.str();
        }

        /// Returns a pair of iterators. Iterating from the begin will filter out interval_types not overlapping @left
        /// and @right
        /// Uses the default end().
        ///
        /// @param left The lower bound of the interval. No intervals strickly to left of this are reached through
        /// iterator
        /// @paran right The upper bound of the interval. No intervals strickly to right of this are reached through
        /// iterator
        /// @return A pair of iterators pointing to the begin and end of query range
        iterator_pair getSubRange( const Interval & queryInterval ) const
        {
            const auto iterators = std::make_pair( begin( queryInterval ), end() );
            return iterators;
        }

        /// Returns a pair of iterators. Iterating from the first element will not filter any interval_types.
        ///
        iterator_pair getFullRange() const { return std::make_pair( begin(), end() ); }

        /// @return Return an iterator pointing to leftmost interval_type of tree.
        iterator begin() const
        {
            Interval queryInterval = m_root->getBoundaries();
            iterator beginIter = begin( queryInterval );
            return beginIter;
        }

        /// @return Return a default end iterator. Used regardless of the filter on the tree.
        iterator end() const { return iterator(); }

    private:
        /// Changes all the left/right boundary values of outmost nodes. Doesn't rebalance tree in current
        /// implementation.
        ///
        /// @param newInterval Resets boundaries to this if needed.
        void widenTreeBoundary( const Interval & newInterval )
        {
            Interval currentBound = m_root->getBoundaries();

            if ( newInterval.start() < currentBound.start() )
            {
                m_root->setLeftBoundary( newInterval );
            }

            if ( newInterval.end() > currentBound.end() )
            {
                m_root->setRightBoundary( newInterval );
            }
        }

        /// @param queryInterval Iterator will skip intervals that do not overlap this.
        /// @return An Iterator pointing to the first overlapping interval in the tree.
        iterator begin( const Interval & queryInterval ) const
        {
            return iterator( m_root->leftMostChild( queryInterval ), queryInterval );
        }

    private:
        node_shared_ptr m_root;  ///< All nodes of binary-tree stem from this

    public:
        /// A forward iterator for the IntervalTree class. Controls interation within and across nodes.
        /// Note: This is unlike standard stl iterators because tree data is not naturally stored sequentially.
        /// This is because the data is not well-ordered. The nodes containing the data are well-ordered using a
        /// depth-first symmetric
        //  traversal.
        /// The data contained on each node is only well-ordered for a given query-interval
        //  (conditioned on where the central point of the node is in relation to the query-interval).
        /// So, the iterator contains a query interval which it uses well-order & filter the data.
        ///
        /// Uses default copy constructor etc.
        class iterator : public std::iterator< std::forward_iterator_tag, interval_type >
        {
        public:
            /// Construct placeholder iterator for end().
            iterator() : m_nodePtr( nullptr ), m_nodeIterator(), m_queryInterval( 0, 1 ) {}

            /// Construct iterator from particular node. It will skip to first data that overlaps the query interval
            explicit iterator( node_ptr nodePtr, const Interval & queryInterval )
                : m_nodePtr( nodePtr ),
                  m_nodeIterator( nodePtr->begin( queryInterval ) ),
                  m_queryInterval( queryInterval )
            {
                this->skipToNextOverlap();
            }

            /// No implementation on how to compare iterators with different query-intervals.
            bool operator==( const iterator & other ) const
            {
                return m_nodePtr == other.m_nodePtr and m_nodeIterator == other.m_nodeIterator;
            }

            /// @return Return true if iterator is not equal to this.
            bool operator!=( const iterator & other ) const { return not( *this == other ); }

            /// Incrementing moves iterator to point to next interval that overlaps.
            // Two iterators may point to the same value but increment diffently depending on their initial range.
            /// Iterating between two iterators with different query-intervals is not well-defined as they would have a
            // different ordering of the tree-node data.
            iterator & operator++()
            {
                ++m_nodeIterator;
                skipToNextOverlap();
                return *this;
            }

            iterator operator++( int )
            {
                iterator tmp = *this;
                operator++();
                return tmp;
            }

            reference operator*() const { return m_nodeIterator.operator*(); }
            reference operator*() { return m_nodeIterator.operator*(); }

            /// Uses the node iterator which points to the actual data.
            pointer operator->() { return m_nodeIterator.operator->(); }

            void moveToBeginOfSubrange( const Interval & newInterval )
            {
                //                assert( newInterval.start() >= m_queryInterval.start() );

                m_queryInterval = newInterval;
                this->moveToNodeBegin();
                this->skipToNextOverlap();
            }

        private:
            /// Skips to first node/node::iterator pair that still overlaps.
            void skipToNextOverlap()
            {
                while ( m_nodePtr != nullptr and m_nodeIterator.isEnd( m_queryInterval ) )
                {
                    m_nodePtr = m_nodePtr->next( m_queryInterval );
                    this->moveToNodeBegin();
                }
            }

            /// Moves parameter m_nodeIterator to point to first interval that overlaps the query interval.
            void moveToNodeBegin()
            {
                if ( m_nodePtr != nullptr )
                {
                    m_nodeIterator = m_nodePtr->begin( m_queryInterval );
                }

                else
                {
                    m_nodeIterator = typename IntervalTreeNode::iterator();
                }
            }

        private:
            node_ptr m_nodePtr;  ///< Node currently iterating on
            typename IntervalTreeNode::iterator m_nodeIterator;
            Interval m_queryInterval;
        };

    private:
        class IntervalTreeNode
        {

        public:
            class iterator;

        public:
            /// Construct an IntervalTreeNode. Node + child-node stores intervals strictly in [left, right).
            // Only intervals containing central point ((left + right)/2) are stored at this level. Rest are stored in
            // descendants.
            ///
            /// @param [in] nodeInterval The interval that the node will span.
            // All intervals that are subsets of this are contained on this node or its children/grandchildren
            /// @param [in] nextAncestor A node_ptr to its first ancestor that is after it. E.g. parent for a left-node.
            explicit IntervalTreeNode( const Interval & nodeInterval, node_ptr nextAncestor )
                : m_nodeInterval( nodeInterval ),
                  m_centre( nodeInterval.centre() ),
                  m_leftChildNode( nullptr ),
                  m_rightChildNode( nullptr ),
                  m_nextAncestor( nextAncestor )
            {
            }

            IntervalTreeNode( const IntervalTreeNode & rhs ) = delete;
            IntervalTreeNode & operator=( const IntervalTreeNode & ) = delete;

            /// @return True if node & descendant nodes contain no intervals.
            bool empty() const
            {
                bool leftIntervalsEmpty = ( m_leftChildNode == nullptr ) ? true : m_leftChildNode->empty();
                bool rightIntervalsEmpty = ( m_rightChildNode == nullptr ) ? true : m_rightChildNode->empty();

                return m_allIntervals.empty() and leftIntervalsEmpty and rightIntervalsEmpty;
            }

            /// @param [in] newIntervalPtr The new pointer to add to the tree.
            /// @param [in] newInterval The interval that the new pointer spans.
            void insert( pointer newIntervalPtr, const Interval & newInterval )
            {
                // interval coming from pure insertion.
                if ( newInterval.start() == m_centre and newInterval.end() == m_centre )
                {
                    m_allIntervals.insert( newIntervalPtr );
                }

                // m_centre in [leftIntervalValue, rightIntervalValue)
                else if ( newInterval.contains( m_centre ) )
                {
                    m_allIntervals.insert( newIntervalPtr );
                    m_intervalsLeftSorted.insert( newIntervalPtr );
                    m_intervalsRightSorted.insert( newIntervalPtr );
                    // Sanity check of input
                    assert( m_intervalsLeftSorted.size() == m_intervalsRightSorted.size() );
                }

                else if ( newInterval.start() > m_centre )
                {
                    if ( m_rightChildNode == nullptr )
                    {
                        // Right-children inherit nextAncestor
                        m_rightChildNode = std::make_shared< IntervalTreeNode >(
                            Interval( m_centre, m_nodeInterval.end() ), m_nextAncestor );
                    }

                    m_rightChildNode->insert( newIntervalPtr, newInterval );
                }

                else if ( newInterval.end() <= m_centre )  // Working with half-intervals.
                {
                    if ( m_leftChildNode == nullptr )
                    {
                        // Left-childre have parent == nextAncestor
                        m_leftChildNode =
                            std::make_shared< IntervalTreeNode >( Interval( m_nodeInterval.start(), m_centre ), this );
                    }

                    m_leftChildNode->insert( newIntervalPtr, newInterval );
                }

                else
                {
                    WECALL_LOG( DEBUG, "Unable to insert interval= " << newInterval.toString() << " On "
                                                                      << this->toString() );
                    throw utils::wecall_exception( "Unable to insert interval" );
                }
            }

            /// Returns the descendant that is only descended through left-children and of greated depth whilst still
            /// having potential
            // to contain overlapping intervals.
            node_ptr leftMostChild( const Interval & queryInterval )
            {
                // All intervals in m_leftChildNode are to left of m_centre by construction.
                if ( m_leftChildNode == nullptr or queryInterval.start() >= m_centre )
                {
                    return this;
                }

                else
                {
                    return m_leftChildNode->leftMostChild( queryInterval );
                }
            }

            /// @return The next node in the tree. For parent, left-child and right-child order is:- left-child, parent
            /// then right-child.
            /// Right-child and all its descendants only contain intervals with left-value >= m_centre by construction.
            //  So, can quickly rule out
            /// These options by going passing to **ancestor in this case.
            /// **This is the uniquely defined ancestor with a centre further to the right. E.g. for a left-child node
            /// its parent.
            node_ptr next( const Interval & queryInterval )
            {
                // All intervals in m_rightChildNode are strickly to right of m_centre by construction.
                if ( m_rightChildNode == nullptr or queryInterval.end() <= m_centre )
                {
                    return m_nextAncestor;
                }

                else
                {
                    return m_rightChildNode->leftMostChild( queryInterval );
                }
            }

            /// @return The nodes outer boundaries. All intervals in this node and decendants are subintervals of this.
            Interval getBoundaries() const { return m_nodeInterval; }

            /// @param [in] newInterval Adjust node right-boundary to this
            void setRightBoundary( const Interval & newInterval )
            {
                m_nodeInterval = Interval( m_nodeInterval.start(), newInterval.end() );

                if ( m_rightChildNode != nullptr )
                {
                    m_rightChildNode->setRightBoundary( newInterval );
                }
            }

            /// @param [in] newInterval Adjust node left-boundary to this.
            void setLeftBoundary( const Interval & newInterval )
            {
                m_nodeInterval = Interval( newInterval.start(), m_nodeInterval.end() );

                if ( m_leftChildNode != nullptr )
                {
                    m_leftChildNode->setLeftBoundary( newInterval );
                }
            }

            /// @return String that uniquely identifies the node within the tree.
            std::string toString() const
            {
                std::stringstream sstrThisID;
                sstrThisID << "Node interval:- " << m_nodeInterval.toString() << ", Centre=" << m_centre;
                return sstrThisID.str();
            }

            /// If the centre of the node is contained in our queryInterval then all intervals at this node intersect
            /// the queryInterval.
            // Default to using left sorted intervals
            /// If the centre is to the right of our queryInterval. Iterate ascendingly (by left_interval_value)
            //  until there is no longer intersection.
            /// If the centre is to the (strict) left of our queryInterval. Iterate descendingly (by
            /// right_interval_value)
            // until there is no longer intersecion.
            ///
            /// @param [in] queryInterval The interval which we are testing for intersection
            iterator begin( const Interval & queryInterval ) const
            {
                if ( queryInterval.overlaps( Interval( m_centre, m_centre ) ) )
                {
                    return iterator( m_allIntervals.cbegin(), m_allIntervals.cend(), iterator::iterator_type::ALL );
                }

                else if ( queryInterval.contains( m_centre ) )
                {
                    return iterator( m_intervalsLeftSorted.cbegin(), m_intervalsLeftSorted.cend(),
                                     iterator::iterator_type::ALL );
                }

                else if ( queryInterval.end() <= m_centre )
                {
                    return iterator( m_intervalsLeftSorted.cbegin(), m_intervalsLeftSorted.cend(),
                                     iterator::iterator_type::FROM_LEFT );
                }

                else if ( queryInterval.start() > m_centre )
                {
                    return iterator( m_intervalsRightSorted.cbegin(), m_intervalsRightSorted.cend(),
                                     iterator::iterator_type::FROM_RIGHT );
                }

                else
                {
                    WECALL_LOG( FATAL, "Error finding overlaps for interval " << queryInterval.toString() << " on "
                                                                               << this->toString() );
                    throw utils::wecall_exception( "Error finding overlaps for interval" );
                }
            }

        private:
            Interval
                m_nodeInterval;  ///< All intervals stored at this node and in decendants are sub-intervals of this.

            const int64_t m_centre;  ///< All intervals stored at this node contain this value.

            node_shared_ptr m_leftChildNode;   ///< Contains all sub-intervals of [m_left, m_centre)
            node_shared_ptr m_rightChildNode;  ///< Contains all sub-intervals of (m_centre, m_right)

            node_ptr m_nextAncestor;  ///< Next node when ascending the tree. If on a left-node = parent.
            // If on right-node = parent's m_nextGrandparent. If on the root = nullptr.

            left_sorted_set m_allIntervals;
            left_sorted_set
                m_intervalsLeftSorted;  ///< Intervals stored ascendingly by left-interval value. No pure insertions
            right_sorted_set
                m_intervalsRightSorted;  ///< Intervals stored descendingly by right-interval value. No pure insertions

        public:
            /// Forward iterator to iterate over either the intervals ordered by left interval value or the
            // intervals reverse ordered by right-interval value.
            /// Want copy and asignment operators to follow be the default
            class iterator : public std::iterator< std::forward_iterator_tag, interval_type >
            {

            public:
                enum class iterator_type : char
                {
                    FROM_LEFT = 1,
                    FROM_RIGHT = 2,
                    ALL = 3
                };

            public:
                iterator() : m_setIterator( NULL ), m_setEnd( NULL ), m_iteratorType( iterator_type::ALL ) {}

                explicit iterator( const set_iterator & internalIterator,
                                   const set_iterator & setEnd,
                                   const iterator_type & iteratorType )
                    : m_setIterator( internalIterator ), m_setEnd( setEnd ), m_iteratorType( iteratorType )
                {
                }

                bool isEnd( const Interval & queryInterval ) const
                {
                    if ( m_iteratorType == iterator_type::ALL )
                    {
                        return m_setIterator == m_setEnd;
                    }

                    else if ( m_iteratorType == iterator_type::FROM_LEFT )
                    {
                        return ( m_setIterator == m_setEnd or
                                 queryInterval.end() <= left_interval_value()( *m_setIterator ) );
                    }

                    else if ( m_iteratorType == iterator_type::FROM_RIGHT )
                    {
                        return ( m_setIterator == m_setEnd or
                                 right_interval_value()( *m_setIterator ) <= queryInterval.start() );
                    }

                    else
                    {
                        return true;
                    }
                }

                /// Uses the underlining multiset iterator's comparison. For looping this is done internally
                bool operator==( const iterator & other ) const { return m_setIterator == other.m_setIterator; }

                bool operator!=( const iterator & other ) const { return not( *this == other ); }

                iterator & operator++()
                {
                    ++m_setIterator;
                    return *this;
                }

                iterator operator++( int )
                {
                    iterator tmp = *this;
                    operator++();
                    return tmp;
                }

                /// Node stores shared_ptrs. So, need to dereference twice to obtain a reference to the object from an
                /// iterator
                reference operator*() const { return **m_setIterator; }
                reference operator*() { return **m_setIterator; }

                /// @return A shared_ptr to an interval contained in the node.
                pointer operator->() { return *m_setIterator; }

                /// @return The type of iteration this is performing
                iterator_type type() const { return m_iteratorType; }

            private:
                set_iterator m_setIterator;
                set_iterator m_setEnd;
                iterator_type m_iteratorType;
            };
        };  // End IntervalTreeNode declaration
    };      // End IntervalTree declaration
}
}

#endif
