// All content Copyright (C) 2018 Genomics plc
#ifndef NEEDLEMAN_WUNSCH_HPP
#define NEEDLEMAN_WUNSCH_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/optional/optional.hpp>

#include "utils/referenceSequence.hpp"
#include "variant/type/variant.hpp"

namespace echidna
{
namespace utils
{
    class NWVariant
    {
    private:
        int64_t m_start;
        int64_t m_end;
        utils::BasePairSequence m_alt;

    public:
        NWVariant( int64_t start, int64_t end, utils::BasePairSequence alt )
            : m_start( start ), m_end( end ), m_alt( alt )
        {
        }

        bool operator==( const NWVariant & A ) const
        {
            return ( this->m_start == A.m_start && this->m_end == A.m_end && this->m_alt == A.m_alt );
        }

        friend std::ostream & operator<<( std::ostream & out, const NWVariant & nWVariant );

        std::shared_ptr< variant::Variant > getVariant( std::string contig,
                                                        int64_t base,
                                                        const utils::referenceSequencePtr_t & referenceSequence ) const
        {
            const auto region = caller::Region( contig, base + m_start, base + m_end );
            return std::make_shared< variant::Variant >( referenceSequence, region, m_alt );
        }
    };

    struct NWPenalties
    {
    private:
        const int32_t m_baseMatch;
        const int32_t m_baseMismatch;
        // A gap of size n incurs a penalty of gapConstant + n times gapLinear + a multiple of gapPosition that depends
        // on the position it starts at
        const int32_t m_insertionConstant;
        const int32_t m_insertionLinear;
        const int32_t m_insertionOpen;
        const int32_t m_deletionConstant;
        const int32_t m_deletionLinear;
        const int32_t m_deletionOpen;
        const int32_t m_indelPosition;

    public:
        NWPenalties( int32_t baseMatchPenalty = 0,
                     int32_t baseMismatchPenalty = -1000,
                     int32_t indelConstantPenalty = -2000,
                     int32_t indelLinearPenalty = -300,
                     int32_t indelPositionPenalty = -1 )
            : m_baseMatch( baseMatchPenalty ),
              m_baseMismatch( baseMismatchPenalty ),
              m_insertionConstant( indelConstantPenalty ),
              m_insertionLinear( indelLinearPenalty ),
              m_insertionOpen( indelConstantPenalty + indelLinearPenalty ),
              m_deletionConstant( indelConstantPenalty ),
              m_deletionLinear( indelLinearPenalty ),
              m_deletionOpen( indelConstantPenalty + indelLinearPenalty ),
              m_indelPosition( indelPositionPenalty )
        {
        }

        NWPenalties( int32_t baseMatchPenalty,
                     int32_t baseMismatchPenalty,
                     int32_t insertionConstantPenalty,
                     int32_t insertionLinearPenalty,
                     int32_t deletionConstantPenalty,
                     int32_t deletionLinearPenalty,
                     int32_t indelPositionPenalty )
            : m_baseMatch( baseMatchPenalty ),
              m_baseMismatch( baseMismatchPenalty ),
              m_insertionConstant( insertionConstantPenalty ),
              m_insertionLinear( insertionLinearPenalty ),
              m_insertionOpen( insertionConstantPenalty + insertionLinearPenalty ),
              m_deletionConstant( deletionConstantPenalty ),
              m_deletionLinear( deletionLinearPenalty ),
              m_deletionOpen( deletionConstantPenalty + deletionLinearPenalty ),
              m_indelPosition( indelPositionPenalty )
        {
        }

        bool operator==( const NWPenalties & A ) const
        {
            return ( this->m_baseMatch == A.m_baseMatch && this->m_baseMismatch == A.m_baseMismatch &&
                     this->m_insertionConstant == A.m_insertionConstant &&
                     this->m_insertionLinear == A.m_insertionLinear && this->m_insertionOpen == A.m_insertionOpen &&
                     this->m_deletionConstant == A.m_deletionConstant && this->m_deletionLinear == A.m_deletionLinear &&
                     this->m_deletionOpen == A.m_deletionOpen && this->m_indelPosition == A.m_indelPosition );
        }

        friend std::ostream & operator<<( std::ostream & out, const NWPenalties & nWPenalties );

        int32_t matchFunction( const char first, const char second ) const
        {
            return ( first != constants::gapChar and first == second ) ? m_baseMatch : m_baseMismatch;
        }

        int32_t insertionOpen() const { return m_insertionOpen; }

        int32_t deletionOpen() const { return m_deletionOpen; }

        int32_t insertionFunction( const std::size_t position, const std::size_t insertionLength ) const
        {
            // The position argument should be the location of the start of the insertion in the reference sequence,
            // where the very start is counted as 0
            return 0 == insertionLength ? 0 : m_insertionConstant + m_insertionLinear * insertionLength +
                                                  m_indelPosition * position;
        }

        int32_t deletionFunction( const std::size_t position, const std::size_t deletionLength ) const
        {
            // The position argument should be the location of the start of the deletion in the alt sequence, where the
            // very start is counted as 0
            return 0 == deletionLength ? 0 : m_deletionConstant + m_deletionLinear * deletionLength +
                                                 m_indelPosition * position;
        }

        int32_t extendInsertion() const { return m_insertionLinear; }

        int32_t extendDeletion() const { return m_deletionLinear; }
    };

    class NeedlemanWunsch
    {
        using scoreMatrix_t = boost::numeric::ublas::matrix< int32_t >;

        using bitType_t = int64_t;
        using traceMatrix_t = boost::numeric::ublas::matrix< bitType_t >;

    public:
        NeedlemanWunsch( const utils::BasePairSequence & referenceString,
                         const utils::BasePairSequence & altString,
                         const NWPenalties & nWPenalties )
            : m_referenceString( referenceString ), m_altString( altString ), m_nWPenalties( nWPenalties )
        {
        }

        int32_t getScoreMatrix();

        std::vector< NWVariant > traceBack() const;

    private:
        struct Backtrace
        {
            std::vector< NWVariant > variants;
            std::size_t refIndex;
            std::size_t altIndex;
            bitType_t trace;
            bool previousIsDeletion;
        };
        scoreMatrix_t m_scoreMatrix;
        traceMatrix_t m_traceMatrix;
        const utils::BasePairSequence & m_referenceString;
        const utils::BasePairSequence & m_altString;
        const NWPenalties & m_nWPenalties;

        const bitType_t m_matchBit = 1;
        const bitType_t m_snpBit = 2;
        const bitType_t m_insertSingleBit = 4;
        const bitType_t m_deleteSingleBit = 8;
        const bitType_t m_insertMultiBit = 16;
        const bitType_t m_deleteMultiBit = 32;
    };
}
}
#endif
