// All content Copyright (C) 2018 Genomics plc
#ifndef CIGAR_HPP
#define CIGAR_HPP

#include <cstdlib>
#include <vector>
#include <string>
#include <memory>

#include "alignment/cigarItems.hpp"

namespace wecall
{
namespace alignment
{

    class CigarItemBuilder
    {
    public:
        static std::shared_ptr< CigarItem > roll( cigarFlags flag, int64_t length );
        static std::shared_ptr< CigarItem > roll( const uint32_t & cigarBitFlag );
        static std::shared_ptr< CigarItem > roll( const std::string & cigarItemString );
    };

    class Cigar
    {
    public:
        /// constructor from pysam bam interface
        Cigar( const uint32_t nOperations, const uint32_t * cigarPtr );

        Cigar( const std::string & cigarString );

        int64_t length() const;
        int64_t lengthInRef() const;
        int64_t lengthInSeq() const;
        int64_t lengthInSeqWithoutSoftClipping() const;

        std::string toString() const;
        int64_t lengthBeforeRefStartPos() const { return m_lengthBeforeRefStartPos; }
        int64_t lengthAfterRefEndPos() const { return m_lengthAfterRefEndPos; }

        referencePositions_t getRefPositions( int64_t startPos ) const;

        utils::Interval getInverseInterval( const utils::Interval & input ) const;

        std::vector< cigarItem_t >::const_iterator cbegin() const { return m_cigarItems.cbegin(); }
        std::vector< cigarItem_t >::const_iterator cend() const { return m_cigarItems.cend(); }

        std::pair< int64_t, int64_t > stripSoftClipping();

        std::size_t size() const { return m_cigarItems.size(); }
        const cigarItem_t & front() const { return m_cigarItems.front(); }
        const cigarItem_t & back() const { return m_cigarItems.back(); }

    private:
        std::vector< cigarItem_t > m_cigarItems;

        int64_t computeLengthBeforeRefStartPos() const;
        int64_t computeLengthAfterRefEndPos() const;

        int32_t m_lengthBeforeRefStartPos;
        int32_t m_lengthAfterRefEndPos;
    };
}  // namespace alignment
}  // namespace wecall

#endif  // CIGAR_HPP
