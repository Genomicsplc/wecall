// All content Copyright (C) 2018 Genomics plc
#ifndef READ_HPP
#define READ_HPP

#include "common.hpp"
#include "io/pysam.hpp"
#include "alignment/cigar.hpp"
#include "utils/sequence.hpp"
#include "variant/type/variant.hpp"
#include "variant/type/breakpoint.hpp"

#include <memory>
#include <string>
#include <cstdint>

namespace wecall
{
namespace io
{

    class Read;
    using readPtr_t = std::shared_ptr< Read >;

    /// The Read class represents a sequence read that has been mapped and aligned to a
    /// reference genome.
    class Read
    {
    public:
        class ReadParams
        {
        public:
            ReadParams();
            ReadParams( const bam1_t * bamRecord );

            ReadParams( std::string seq,
                        utils::QualitySequence qual,
                        std::string readGroupID,
                        alignment::Cigar cigar,
                        int32_t tid,
                        int32_t startPos,
                        uint16_t flag,
                        uint8_t mappingQuality,
                        int32_t insertSize,
                        int32_t mateTid,
                        int32_t mateStartPos,
                        const std::string & qname );

            std::string m_sequenceStr;
            utils::QualitySequence m_qualities;
            std::string m_qname;
            std::string m_readGroupID;

            alignment::Cigar m_cigar;

            int32_t m_tid;
            int32_t m_startPos;  // Start position of read in reference
            int32_t m_endPos;    // End position of read in reference
            uint16_t m_flag;
            uint8_t m_mappingQuality;
            int32_t m_insertSize;
            int32_t m_mateTid;
            int32_t m_mateStartPos;  // Start position of mate, for paired-end reads.
        };

        Read( const ReadParams & params, utils::referenceSequencePtr_t refSequence );

        /// Construct a Read from a Samtools BAM record
        ///
        /// @param a Samtools BAM record
        Read( const bam1_t * bamRecord, utils::referenceSequencePtr_t refSequence );

        /// Construct a Read from data,
        ///
        /// @param seq The read seqeunce
        /// @param qual The read quality scores, as PHRED values
        Read( utils::BasePairSequence seq,
              utils::QualitySequence qual,
              std::string readGroupID,
              alignment::Cigar cigar,
              int32_t tid,
              int32_t startPos,
              uint16_t flag,
              uint8_t mappingQuality,
              int32_t insertSize,
              int32_t mateTid,
              int32_t mateStartPos,
              utils::referenceSequencePtr_t refSequence,
              std::string qname = std::string() );

        /// Move constructor
        ///
        /// @param rhs Read to be moved
        Read( Read && rhs ) = default;

        /// Move assignment operator
        ///
        /// @param rhs Read to be moved
        /// @return reference to this read
        Read & operator=( Read && rhs ) = default;

        /// Return a string representation of this read
        std::string toString() const;

        /// Return the sequence of As, Cs, Ts and Gs
        const utils::BasePairSequence & sequence() const
        {
            if ( m_sequence.size() == 0 and m_isReference )
            {
                const_cast< utils::BasePairSequence & >( m_sequence ) = makeRefSequence();
            }
            return m_sequence;
        }

        /// Return the unique name (usually refered to as qname) of the read.
        const std::string & getQName() const { return m_qname; }

        /// Return the sequence of quality scores (one per base) as PHRED scores.
        const utils::QualitySequence & getQualities() const { return m_qualities; }
        const alignment::Cigar & cigar() const { return m_cigar; }

        // Allows to modify qualities in-place and reuse the rest of this data structure
        // if you change qualities tag this read appropriately; presently we only
        // modify qualities through error corrector so don't support tags but rather
        // call setErrorCorrectorModifiedQualities() through which we also communicate
        // error-corrector specific information about where error corrector kicked in
        utils::QualitySequence & qualities() { return m_qualities; }

        /// Return the start position of the read
        int64_t getStartPos() const { return m_startPos; }

        caller::Region getRegion() const
        {
            return caller::Region( m_refSequence->contig(), m_startPos, m_alignedEndPos );
        }

        // TODO(ES): Add abilty to retrieve the mate's contig.
        bool isMateOnSameContig() const { return m_tid == m_mateTid; }

        // TODO(ES): Use insert size to make this exact.
        utils::Interval getMateIntervalInRef() const
        {
            return utils::Interval( m_mateStartPos, m_mateStartPos + getAlignedLength() );
        }

        utils::Interval getMaximalReadInterval() const;

        /// Return the end position of the read (eqivalent to startPos + readLength)
        int64_t getEndPos() const { return m_endPos; }

        /// Return the flag from the original BAM record
        int64_t getFlag() const { return m_flag; }

        /// Return the aligned end position of the read. This is not always equal to the value
        /// returned by getEndPos as it takes into account alignment gaps i.e. insertions and
        /// deletions
        int64_t getAlignedEndPos() const { return m_alignedEndPos; }

        /// Return the amount of read that is before the start position. The start position is taken from the reference.
        /// If the cigar starts with insertions then these will be before the (aligned) start position.
        int64_t getLengthBeforeAlignedStartPos() const { return m_cigar.lengthBeforeRefStartPos(); }

        /// Return the amount of read that is after the aligned end position. The aligned end position refers to the
        /// reference.
        /// If the cigar ends with insertions then these will be after the aligned end position.
        int64_t getLengthAfterAlignedEndPos() const;

        /// Return the length of the read
        int64_t getLength() const { return m_endPos - m_startPos; }

        /// Return the aligned length of the read on the reference sequence, taking into account
        /// insertions and deletions
        int64_t getAlignedLength() const { return m_alignedEndPos - m_startPos; }

        /// Return the length of the original read-pair fragment, as inferred by the read-mapper from
        /// the start positions of read1 and read2. This can be negative.
        int64_t getInsertSize() const { return m_insertSize; }

        /// Return the start position of the mate of this read
        int64_t getMateStartPos() const { return m_mateStartPos; }

        /// Return mapping quality as a PHRED score
        int64_t getMappingQuality() const { return m_mappingQuality; }

        /// Return the ID of this read as listed in the @RG tag.
        std::string getReadGroupID() const { return m_readGroupID; }

        /// Is this read in the reverse direction (the reference sequence is always forward, by
        /// convention).
        bool isReverse() const;

        /// Is this read one of a pair
        bool isPaired() const;

        /// Has this read been marked as a duplicate e.g. by Samtools or PICARD
        bool isDuplicate() const;

        /// Is this secondary alignment
        bool isSecondary() const;

        /// Return true if the read is not mapped. Unmapped reads with mapped mates are placed next
        /// to their mates in the BAM format, whilst pairs of unmapped reads go at the end of the file
        bool isUnMapped() const;

        /// Is the mate of this read unmapped
        bool isMateUnMapped() const;

        /// Is the mate of this read in the reverse direction
        bool isMateReverse() const;

        /// Is this read the first read of the pair
        bool isReadOne() const;

        /// Return true if each segment is properly aligned to according to the aligner.
        bool isProperPair() const;

        // matches the reference genome.
        bool isReference() const { return m_isReference; }

        /// Trim overlapping part of forward read, in pairs where the read length is greater than the insert size
        /// N.B Insert size is from start of forward read to end of reverse read, i.e. fragment size. This is done to
        /// remove duplicate information, which gives systematic errors when pcr errors have occured in library prep.
        void trimOverlap();

        /// Trim the end of any read where the insert size is < read length. If these have not been
        /// already filtered out then they need trimming, as adapter contamination will cause a
        /// high FP rate otherwise.
        void trimReadOfShortFragment();

        /// For each character in the sequence return a the corresponding position in the reference.
        /// Uses emptyPos of -1 for inserted sequence.
        alignment::referencePositions_t getReferencePositions() const;

        // Take an interval aligned to the reference and return corresponding interval in read space.
        utils::Interval getIntervalInRead( const utils::Interval & refInterval ) const;

        std::vector< variant::varPtr_t > getVariants() const;

        std::vector< variant::breakpointPtr_t > getBreakpoints() const;

    private:
        Read( const Read & rhs ) = delete;

        utils::BasePairSequence makeRefSequence() const;
        std::pair< utils::BasePairSequence::const_iterator, utils::BasePairSequence::const_iterator >
        getRefSequenceRange() const;

    private:
        const utils::BasePairSequence m_sequence;
        utils::QualitySequence m_qualities;
        std::string m_qname;
        std::string m_readGroupID;

        alignment::Cigar m_cigar;

        const int32_t m_tid;
        const int32_t m_startPos;       // Start position of read in reference
        const int32_t m_endPos;         // End position of read in reference
        const int32_t m_alignedEndPos;  // Actual end of read in reference, when alignment is taken into account
        const int32_t m_insertSize;
        const int32_t m_mateTid;
        const int32_t m_mateStartPos;  // Start position of mate, for paired-end reads.

        utils::referenceSequencePtr_t m_refSequence;
        const uint16_t m_flag;
        const uint8_t m_mappingQuality;
        uint8_t m_isReference = false;
    };
}
}

#endif
