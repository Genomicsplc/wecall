// All content Copyright (C) 2018 Genomics plc
#ifndef CIGAR_ITEMS_HPP
#define CIGAR_ITEMS_HPP

#include <vector>
#include <string>
#include <memory>
#include "common.hpp"
#include "variant/variantGenerationData.hpp"
#include "variant/type/variant.hpp"

namespace wecall
{
namespace variant
{
    class Variant;
    using varPtr_t = std::shared_ptr< Variant >;
}

namespace alignment
{
    /// Enum of CIGAR string flags. The CIGAR string stores pairs of
    /// integers, (flag, length) describing the alignment of a read to
    /// the reference sequence. This is a utility for mapping the flag
    /// values to descriptive names.
    enum class cigarFlags
    {
        MATCH = 0,        // Match
        INSERTION = 1,    // Insertion
        DELETION = 2,     // Deletion
        SKIP = 3,         // Skipped region from reference
        SOFT_CLIP = 4,    // Soft clipping. Sequence is present in read
        HARD_CLIP = 5,    // Hard clipping. Sequence is not present in read
        PADDING = 6,      // Padding. Used for padded alignment
        SEQ_MATCH = 7,    // Sequence match.
        SEQ_MISMATCH = 8  // Sequence mismatch.
    };

    const int emptyPos = -1;  ///< Denotes lack of alignment with reference
    const int noPos = -1;     ///< Denotes no start and/or end in collection of reads (collection is empty).

    using referencePositions_t = std::vector< int64_t >;

    //-----------------------------------------------------------------------------------------

    struct Offsets
    {
        Offsets( int64_t readOffset, int64_t refOffset ) : read( readOffset ), ref( refOffset ) {}
        int64_t read;
        int64_t ref;
    };

    using offsetsPtr_t = std::shared_ptr< Offsets >;

    //-----------------------------------------------------------------------------------------

    class CigarItem
    {
    public:
        CigarItem( int length ) : m_length( length ) {}

        virtual ~CigarItem() {}

        int64_t length() const { return m_length; }

        virtual int64_t lengthInRef() const = 0;
        virtual int64_t lengthInSeq() const = 0;
        void moveOffsets( offsetsPtr_t offsets ) const;
        virtual referencePositions_t getRelativeRefPositions( const int64_t & startRefPos ) const = 0;

        virtual variant::variantSet_t getVariants( variant::variantGenerationDataPtr_t variantGenerationData,
                                                   offsetsPtr_t offsets ) const
        {
            return {};
        }

        std::string toString() const { return std::to_string( m_length ) + this->operationSymbol(); }

        virtual bool isSoftClipped() const { return false; }
        virtual bool isHardClipped() const { return false; }

    protected:
        virtual char operationSymbol() const = 0;

        int64_t m_length;
    };

    using cigarItem_t = std::shared_ptr< CigarItem >;

    //-----------------------------------------------------------------------------------------

    class CigarMatch : public CigarItem
    {
    public:
        CigarMatch( int length ) : CigarItem( length ) {}

        virtual ~CigarMatch() {}

        int64_t lengthInRef() const override { return this->length(); }
        int64_t lengthInSeq() const override { return this->length(); }
        referencePositions_t getRelativeRefPositions( const int64_t & startRefPos ) const override;

        variant::variantSet_t getVariants( variant::variantGenerationDataPtr_t variantGenerationData,
                                           offsetsPtr_t offsets ) const override;

    private:
        char operationSymbol() const override { return 'M'; }
    };

    //-----------------------------------------------------------------------------------------

    class CigarInsertion : public CigarItem
    {
    public:
        CigarInsertion( int length ) : CigarItem( length ) {}

        virtual ~CigarInsertion() {}

        int64_t lengthInRef() const override { return 0; }
        int64_t lengthInSeq() const override { return this->length(); }
        referencePositions_t getRelativeRefPositions( const int64_t & startRefPos ) const override;

        variant::variantSet_t getVariants( variant::variantGenerationDataPtr_t varGenData,
                                           offsetsPtr_t offsets ) const override;

    private:
        char operationSymbol() const override { return 'I'; }
    };

    //-----------------------------------------------------------------------------------------

    class CigarDeletion : public CigarItem
    {
    public:
        CigarDeletion( int length ) : CigarItem( length ) {}

        virtual ~CigarDeletion() {}

        int64_t lengthInRef() const override { return this->length(); }
        int64_t lengthInSeq() const override { return 0; }
        referencePositions_t getRelativeRefPositions( const int64_t & startRefPos ) const override { return {}; }

        variant::variantSet_t getVariants( variant::variantGenerationDataPtr_t varGenData,
                                           offsetsPtr_t offsets ) const override;

    private:
        char operationSymbol() const override { return 'D'; }
    };

    //-----------------------------------------------------------------------------------------

    class CigarSkip : public CigarItem
    {
    public:
        CigarSkip( int length ) : CigarItem( length ) {}

        virtual ~CigarSkip() {}

        int64_t lengthInRef() const override { return this->length(); }
        int64_t lengthInSeq() const override { return 0; }
        referencePositions_t getRelativeRefPositions( const int64_t & startRefPos ) const override { return {}; }

    private:
        char operationSymbol() const override { return 'N'; }
    };

    //-----------------------------------------------------------------------------------------

    class CigarSoftClip : public CigarItem
    {
    public:
        CigarSoftClip( int length ) : CigarItem( length ) {}

        virtual ~CigarSoftClip() {}

        int64_t lengthInRef() const override { return 0; }
        int64_t lengthInSeq() const override { return this->length(); }
        referencePositions_t getRelativeRefPositions( const int64_t & startRefPos ) const override;
        virtual bool isSoftClipped() const override { return true; }

    private:
        char operationSymbol() const override { return 'S'; }
    };

    //-----------------------------------------------------------------------------------------

    class CigarHardClip : public CigarItem
    {
    public:
        CigarHardClip( int length ) : CigarItem( length ) {}

        virtual ~CigarHardClip() {}

        int64_t lengthInRef() const override { return 0; }
        int64_t lengthInSeq() const override { return 0; }
        referencePositions_t getRelativeRefPositions( const int64_t & startRefPos ) const override { return {}; }
        virtual bool isHardClipped() const override { return true; }

    private:
        char operationSymbol() const override { return 'H'; }
    };

    //-----------------------------------------------------------------------------------------

    class CigarPadding : public CigarItem
    {
    public:
        CigarPadding( int length ) : CigarItem( length ) {}

        virtual ~CigarPadding() {}

        int64_t lengthInRef() const override { return 0; }
        int64_t lengthInSeq() const override { return 0; }
        referencePositions_t getRelativeRefPositions( const int64_t & startRefPos ) const override { return {}; }

    private:
        char operationSymbol() const override { return 'P'; }
    };

    //-----------------------------------------------------------------------------------------

    class CigarSequenceMatch : public CigarMatch
    {
    public:
        CigarSequenceMatch( int length ) : CigarMatch( length ) {}
        virtual ~CigarSequenceMatch() {}

    private:
        char operationSymbol() const override { return '='; }
    };

    //-----------------------------------------------------------------------------------------

    class CigarSequenceMismatch : public CigarMatch
    {
    public:
        CigarSequenceMismatch( int length ) : CigarMatch( length ) {}
        virtual ~CigarSequenceMismatch() {}

    private:
        char operationSymbol() const override { return 'X'; }
    };

}  // namespace alignment
}  // namespace wecall

#endif  // CIGAR_ITEMS_HPP
