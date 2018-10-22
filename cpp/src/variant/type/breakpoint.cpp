// All content Copyright (C) 2018 Genomics plc
#include "common.hpp"
#include "variant/type/breakpoint.hpp"
#include "caller/region.hpp"
#include <sstream>
#include <variant/haplotype.hpp>

namespace wecall
{
namespace variant
{
    Breakpoint::Breakpoint( const std::string & contig,
                            const int64_t pos,
                            const bool isStartBreakpoint,
                            const utils::BasePairSequence & sequence )
        : m_contig( contig ), m_isStartBreakpoint( isStartBreakpoint ), m_pos( pos ), m_sequence( sequence )
    {
    }

    bool Breakpoint::operator==( const Breakpoint & other ) const
    {
        return m_contig == other.m_contig and m_isStartBreakpoint == other.m_isStartBreakpoint and
               m_pos == other.m_pos and m_sequence == other.m_sequence;
    }
    void Breakpoint::addLocalVariants( const variantSet_t & localVariants )
    {
        m_localVariants.insert( localVariants.cbegin(), localVariants.cend() );
    }

    void Breakpoint::addLocalVariant( const varPtr_t & var ) { m_localVariants.insert( var ); }

    void Breakpoint::addQueryRegion( const caller::Region & region ) { m_mateRegions.insert( region ); }

    std::string Breakpoint::toString() const
    {
        std::stringstream sstr;
        sstr << ( isStartBreakpoint() ? "Srt" : "End" );
        sstr << "Breakpoint(";
        // Create a nice representation here.
        sstr << contig() << ":" << m_pos;
        sstr << " " << sequence();
        sstr << ")";
        sstr << "\n\tLocal variants:\t";
        for ( const auto & var : this->getLocalVariants() )
        {
            sstr << *var << ", ";
        }
        sstr << "\n\tMateRegions:\t";
        for ( const auto & region : m_mateRegions )
        {
            sstr << region << ", ";
        }
        return sstr.str();
    }

    bool breakpointPtrComp::operator()( const breakpointPtr_t & x, const breakpointPtr_t & y ) const
    {
        if ( x->contig() != y->contig() )
        {
            return x->contig() < y->contig();
        }
        if ( x->pos() != y->pos() )
        {
            return x->pos() < y->pos();
        }
        if ( x->isStartBreakpoint() != y->isStartBreakpoint() )
        {
            return y->isStartBreakpoint();
        }
        if ( x->sequence().size() != y->sequence().size() )
        {
            return x->sequence().size() < y->sequence().size();
        }
        return x->sequence() < y->sequence();
    }

    bool BreakpointLocus::hasMinSupport() const
    {
        for ( const auto & pair : m_breakpoints )
        {
            if ( pair.second >= m_minSupport )
            {
                return true;
            }
        }
        return false;
    }

    void BreakpointLocus::add( breakpointPtr_t newBp )
    {
        WECALL_ASSERT(
            m_contig == newBp->contig() and m_pos == newBp->pos() and m_isStartLocus == newBp->isStartBreakpoint(),
            "" );

        if ( this->isStartLocus() and newBp->sequence().size() >= m_minStartSize )
        {
            ++m_breakpoints[newBp->sequence().substr( 0, m_minStartSize )];
        }
        else if ( newBp->sequence().size() >= m_minStartSize )
        {
            ++m_breakpoints[newBp->sequence().substr( newBp->sequence().size() - m_minStartSize, m_minStartSize )];
        }

        const auto localVariants = newBp->getLocalVariants();
        m_localVariants.insert( localVariants.cbegin(), localVariants.cend() );

        for ( const auto & region : newBp->mateRegions() )
        {
            this->addMateRegion( region );
        }
    }

    void BreakpointLocus::addMateRegion( const caller::Region & mateRegion, const std::size_t padding )
    {
        if ( ( this->isStartLocus() and mateRegion.end() <= pos() ) or
             ( not this->isStartLocus() and mateRegion.start() >= pos() ) )
        {
            return;
        }

        const caller::Region paddedRegion( mateRegion.contig(),
                                           ( this->isStartLocus() ? mateRegion.start() - padding : mateRegion.start() ),
                                           ( this->isStartLocus() ? mateRegion.end() : mateRegion.end() + padding ) );

        m_mateRegions.insert( paddedRegion );
    }

    std::string BreakpointLocus::toString() const
    {
        std::stringstream sstr;
        if ( this->isStartLocus() )
        {
            sstr << "Start";
        }
        else
        {
            sstr << "End";
        }
        sstr << "BreakpointLocus(" << this->contig() << ":" << this->pos() << "\t";
        if ( not m_mateRegions.empty() )
        {
            sstr << "MateRegions: " << m_mateRegions.getSpan();
        }
        sstr << ")";
        return sstr.str();
    }

    bool BreakpointLocusCluster::isRelated( const breakpointLocusPtr_t & breakpointLocusPtr ) const
    {
        const auto bpRegion = caller::Region( breakpointLocusPtr->contig(), breakpointLocusPtr->pos(),
                                              breakpointLocusPtr->pos() ).getPadded( m_paddingDistance );
        bool related = false;
        if ( m_regions.overlaps( bpRegion ) )
        {
            related = true;
        }
        for ( const auto & bpMate : breakpointLocusPtr->mateRegions() )
        {
            if ( m_regions.overlaps( bpMate ) )
            {
                related = true;
            }
        }

        return related;
    }

    void BreakpointLocusCluster::push( const breakpointLocusPtr_t & breakpointLocusPtr )
    {
        m_breakpoints.push_back( breakpointLocusPtr );
        for ( const auto & region : breakpointLocusPtr->mateRegions() )
        {
            m_regions.insert( region );
        }
        m_regions.insert( caller::Region( breakpointLocusPtr->contig(), breakpointLocusPtr->pos(),
                                          breakpointLocusPtr->pos() ).getPadded( m_paddingDistance ) );
        m_regions.fill( m_paddingDistance );
    }

    bool breakpointLocusPtrComp::operator()( const breakpointLocusPtr_t & x, const breakpointLocusPtr_t & y ) const
    {
        if ( x->contig() != y->contig() )
        {
            return x->contig() < y->contig();
        }
        else if ( x->pos() != y->pos() )
        {
            return x->pos() < y->pos();
        }
        else if ( x->isStartLocus() != y->isStartLocus() )
        {
            return not x->isStartLocus();
        }
        else
        {
            return false;
        }
    }
}
}
