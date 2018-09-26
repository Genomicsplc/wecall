// All content Copyright (C) 2018 Genomics plc
#ifndef VARIANT_HPP
#define VARIANT_HPP

#include <cstdint>
#include <cstdlib>
#include <memory>
#include <ostream>
#include <tuple>
#include <string>
#include <vector>
#include <set>
#include <boost/optional.hpp>
#include "stats/functions.hpp"
#include "utils/referenceSequence.hpp"

#include "caller/region.hpp"
#include "utils/sequence.hpp"

namespace echidna
{
namespace io
{

    class Read;
    using readPtr_t = std::shared_ptr< Read >;
}
}

namespace echidna
{
namespace variant
{

    //-------------------------------------------------------------------------------------------------
    // The Variant class encapsulates information for a specific genetic mutation, e.g.
    // a single SNP, indel, structural variant or other variant type. Here we define a pure abstract
    // base class that all variant classes must implement.
    //-------------------------------------------------------------------------------------------------

    class Variant;

    using varPtr_t = std::shared_ptr< Variant >;

    /// Comparison function used to sort shared pointers to variant objects. Use this with sort
    /// to sort a vector/set/whatever of variant pointers.
    struct varPtrComp
    {
        bool operator()( const varPtr_t & x, const varPtr_t & y ) const;
    };

    using variantSet_t = std::set< varPtr_t, varPtrComp >;

    class Variant
    {
    public:
        typedef utils::Interval Interval;
        typedef utils::ReferenceSequence ReferenceSequence;
        typedef utils::BasePairSequence BasePairSequence;
        typedef caller::Region Region;

        Variant( utils::referenceSequencePtr_t referenceSequencePtr,
                 Region region,
                 BasePairSequence added,
                 bool isFullyLeftAligned = false );

        Variant operator-( const Variant & other ) const;
        Variant operator+( const Variant & other ) const;
        bool operator==( const Variant & other ) const;
        friend std::ostream & operator<<( std::ostream & out, const Variant & var ) { return out << var.toString(); }

        std::string contig() const { return m_region.contig(); }
        int64_t start() const { return m_region.start(); }
        int64_t end() const { return m_region.end(); }

        Interval interval() const { return region().interval(); }
        Region region() const { return m_region; }

        std::string toString() const;

        int64_t sequenceLength() const { return static_cast< int64_t >( m_altSequence.size() ); }
        int64_t sequenceLengthInRef() const { return static_cast< int64_t >( m_region.size() ); }
        int64_t sequenceLengthChange() const { return sequenceLength() - sequenceLengthInRef(); }

        bool empty() const { return sequenceLengthInRef() == 0 and sequenceLength() == 0; }

        bool isPureIndel() const { return this->sequenceLength() == 0 or this->sequenceLengthInRef() == 0; }
        bool isIndel() const { return 0 != this->sequenceLengthChange(); }
        bool isDeletion() const { return this->sequenceLengthInRef() > this->sequenceLength(); }
        bool isInsertion() const { return this->sequenceLengthInRef() < this->sequenceLength(); }
        bool isSNP() const { return this->sequenceLength() == 1 and this->sequenceLengthInRef() == 1; }
        bool isMNP() const
        {
            return this->sequenceLength() == this->sequenceLengthInRef() and this->sequenceLength() > 1;
        }
        bool isLargeVariant( const int largeVariantThreshold ) const
        {
            return this->sequenceLengthInRef() >= largeVariantThreshold;
        }

        const BasePairSequence & sequence() const { return m_altSequence; }
        ReferenceSequence refSequence() const { return m_refSequence->subseq( m_region ); }

        bool overlaps( Region region ) const { return m_region.overlaps( region ); }
        bool overlaps( varPtr_t other ) const { return this->overlaps( other->region() ); };

        varPtr_t getLeftAligned( const int64_t minPos ) const;
        varPtr_t getRightAligned( const int64_t maxPos ) const;
        caller::Region getLeftAlignedRegion( const boost::optional< int64_t > minPos = boost::none ) const;
        caller::Region getRightAlignedRegion( const boost::optional< int64_t > maxPos = boost::none ) const;
        bool isFullyLeftAligned() const { return m_isFullyLeftAligned; }

        bool removable( const varPtr_t & other ) const;
        varPtr_t remove( const varPtr_t & other ) const;
        bool joinable( const varPtr_t & other ) const;
        varPtr_t join( const varPtr_t & other ) const { return ( *this + *other ).getTrimmed(); }

        varPtr_t getTrimmed() const;
        std::vector< varPtr_t > split() const;

        caller::SetRegions getStartEndRegions( const boost::optional< int64_t > minPos = boost::none,
                                               const boost::optional< int64_t > maxPos = boost::none ) const;

        /// This function is defined purely for convenience, as it's easier to remember and
        /// more intuitive
        int64_t zeroIndexedVcfPosition() const { return this->isPureIndel() ? m_region.start() - 1 : m_region.start(); }

        bool neverFilter() const { return m_neverFilter; }
        void disableFiltering() { m_neverFilter = true; }

        bool isGenotypingVariant() const { return m_isGenotypingVariant; }
        void setGenotypingVariant() { m_isGenotypingVariant = true; }

        bool isFromBreakpoint() const { return m_fromBreakpoint; }
        void setFromBreakpoint() { m_fromBreakpoint = true; }

        double prior() const { return m_prior; };
        phred_t phredScaledPrior() const { return stats::toPhredQ( this->prior() ); }
        void prior( const double prior ) { m_prior = std::max( prior, constants::minVariantPrior ); }

        std::vector< io::readPtr_t > getReads() const { return m_reads; }
        void addRead( io::readPtr_t readPtr ) { m_reads.push_back( readPtr ); }

    private:
        utils::referenceSequencePtr_t m_refSequence;
        const Region m_region;
        const BasePairSequence m_altSequence;

        const bool m_isFullyLeftAligned;

        double m_prior;

        bool m_neverFilter;
        bool m_isGenotypingVariant;
        bool m_fromBreakpoint;

        std::vector< io::readPtr_t > m_reads;
    };

    void setDefaultPriors( const std::vector< varPtr_t > & variants );
}
}

#endif
