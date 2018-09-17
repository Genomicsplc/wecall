// All content Copyright (C) 2018 Genomics plc
#ifndef BREAKPOINT_HPP
#define BREAKPOINT_HPP
#include <boost/optional.hpp>
#include <ostream>
#include <set>
#include <map>

#include "common.hpp"
#include "variant/type/variant.hpp"
#include "utils/sequence.hpp"
#include "caller/region.hpp"
#include "utils/sequence.hpp"
#include "utils/logging.hpp"

namespace echidna
{
namespace variant
{
    class Breakpoint
    {
    public:
        Breakpoint( const std::string & contig,
                    const int64_t pos,
                    const bool isStartBreakpoint,
                    const utils::BasePairSequence & sequence );

        friend std::ostream & operator<<( std::ostream & out, const Breakpoint & var ) { return out << var.toString(); }
        bool operator==( const Breakpoint & other ) const;

        std::string contig() const { return m_contig; }
        int64_t pos() const { return m_pos; }
        bool isStartBreakpoint() const { return m_isStartBreakpoint; }

        const utils::BasePairSequence & sequence() const { return m_sequence; }

        std::string toString() const;

        void addQueryRegion( const caller::Region & region );

        const caller::SetRegions & mateRegions() const { return m_mateRegions; }

        void addLocalVariants( const variantSet_t & localVariants );
        void addLocalVariant( const varPtr_t & localVariants );
        variantSet_t getLocalVariants() const { return m_localVariants; }

    private:
        // To mirror caller::Region object except with optional start or end.
        const std::string m_contig;
        // Convert to enum
        const bool m_isStartBreakpoint;
        const int64_t m_pos;

        const utils::BasePairSequence m_sequence;

        caller::SetRegions m_mateRegions;
        variantSet_t m_localVariants;
    };

    typedef std::shared_ptr< Breakpoint > breakpointPtr_t;

    struct breakpointPtrComp
    {
        bool operator()( const breakpointPtr_t & x, const breakpointPtr_t & y ) const;
    };

    class BreakpointLocus
    {
    public:
        BreakpointLocus( const std::string & contig, const int64_t pos, const bool isStartBreakpoint )
            : m_contig( contig ), m_pos( pos ), m_isStartLocus( isStartBreakpoint ){};

        friend std::ostream & operator<<( std::ostream & out, const BreakpointLocus & bpLocus )
        {
            return out << bpLocus.toString();
        }

        std::string contig() const { return m_contig; }
        int64_t pos() const { return m_pos; }
        bool isStartLocus() const { return m_isStartLocus; }

        void add( breakpointPtr_t newBp );
        bool hasMinSupport() const;
        std::string toString() const;

        void addMateRegion( const caller::Region & mateRegion, const std::size_t padding = 150 );
        const caller::SetRegions & mateRegions() const { return m_mateRegions; }

        variantSet_t getLocalVariants() const { return m_localVariants; }

    private:
        const std::string m_contig;
        const int64_t m_pos;
        const bool m_isStartLocus;

        std::map< utils::BasePairSequence, std::size_t > m_breakpoints;
        caller::SetRegions m_mateRegions;
        variantSet_t m_localVariants;

        const std::size_t m_minStartSize = 3;
        const std::size_t m_minSupport = 2;
    };

    typedef std::shared_ptr< BreakpointLocus > breakpointLocusPtr_t;

    struct breakpointLocusPtrComp
    {
        bool operator()( const breakpointLocusPtr_t & x, const breakpointLocusPtr_t & y ) const;
    };

    typedef std::set< breakpointLocusPtr_t, breakpointLocusPtrComp > breakpointLocusSet_t;

    class BreakpointLocusCluster
    {
    public:
        BreakpointLocusCluster( const std::size_t paddingDistance ) : m_paddingDistance( paddingDistance ) {}

        bool isRelated( const breakpointLocusPtr_t & breakpointLocusPtr ) const;
        void push( const breakpointLocusPtr_t & breakpointLocusPtr );

        const std::vector< breakpointLocusPtr_t > & get() const { return m_breakpoints; }

    private:
        const std::size_t m_paddingDistance;
        std::vector< breakpointLocusPtr_t > m_breakpoints;
        caller::SetRegions m_regions;
    };
}
}

#endif
