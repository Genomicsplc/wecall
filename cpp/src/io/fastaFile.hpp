// All content Copyright (C) 2018 Genomics plc
#ifndef FASTA_FILE_HPP
#define FASTA_FILE_HPP

#include "common.hpp"
#include "utils/referenceSequence.hpp"
#include "utils/logging.hpp"
#include "utils/timer.hpp"
#include "utils/sequence.hpp"
#include "caller/region.hpp"

#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <memory>
#include <map>
#include <set>

namespace echidna
{
namespace io
{
    using std::ifstream;

    /// A simple struct that stores information from a single line of the
    /// FASTA index file. All data is public.
    struct IndexTuple
    {
    public:
        /// Constructor
        ///
        /// @param seqName The name of the chromosome/contig
        /// @param seqLength The length of the sequence
        /// @param start The starting position, in the file, of the sequence
        /// @param lineLength The length of a single line in the FASTA file
        /// @param fullLineLength The length of a line including newlines etc.
        IndexTuple( const std::string & seqName,
                    const int64_t seqLength,
                    const int64_t start,
                    const int64_t lineLength,
                    const int64_t fullLineLength )
            : m_seqName( seqName ),
              m_startPos( start ),
              m_seqLength( seqLength ),
              m_lineLength( lineLength ),
              m_fullLineLength( fullLineLength )
        {
        }

        const std::string m_seqName;
        const int64_t m_startPos;
        const int64_t m_seqLength;
        const int64_t m_lineLength;
        const int64_t m_fullLineLength;
    };

    using itPtr_t = std::shared_ptr< IndexTuple >;
    /// Storage container for .fai index information.
    typedef std::map< std::string, itPtr_t > contigMap_t;

    std::string fastaIndexFileName( std::string fastaFileName );

    /// Deals with reading the .fai index file and working out the offsets
    /// of the various contigs.
    class FastaIndex
    {
    public:
        /// Constructor
        ///
        /// @param fileName The name of the FASTA file
        explicit FastaIndex( const std::string & fileName );

        /// Destructor
        ~FastaIndex() {}

        /// Retrieve a pointer to the IndexTuple instance that represents this
        /// reference contig.
        ///
        /// @param refName The name of the contig to retrieve
        /// @return A pointer to the IndexTuple
        const itPtr_t getIndexTuple( const std::string & refName ) const;

        /// Retrieve a vecto of the names of all the contigs/chromosomes stored in the FASTA file
        ///
        /// @return A vector of strings of the contig names
        std::vector< std::string > standardContigs() const { return m_standardContigs; }

        /// Retrieves start position of contig in ref file. Will throw out_of_range exception if contig is not
        /// present.
        ///
        int64_t contigStart( const std::string & contig ) const { return m_contigMap.at( contig )->m_startPos; }

        int64_t getContigLength( const std::string & refName ) const
        {
            return this->getIndexTuple( refName )->m_seqLength;
        }

        std::map< std::string, utils::Interval > contigs() const;

    private:
        /// Read the index file and extract information about the size and location in the FASTA
        /// file of each contig.
        ///
        /// @param theFile The index file
        void parseContigInfoFromFile( ifstream & theFile );

        contigMap_t m_contigMap;

        std::vector< std::string > m_standardContigs;  ///< _ordered_ set of standard contigs

        utils::timerPtr_t m_timer;
    };

    /// A representation of a FASTA file. Maintains an in-memory cache of a large chunk
    /// of sequence, which makes sequence retrieval much faster than reading from the file each time.
    class FastaFile
    {
    public:
        /// Constructor.
        ///
        /// @param fileName The name of the FASTA file to read
        explicit FastaFile( const std::string & fileName );

        /// Return a string of the reference sequence. If ends positions are outside of the reference, pad with "N"
        /// characters
        utils::ReferenceSequence getSequence( const caller::Region & region ) const;

        /// Set the internal cache to store the sequence of region.
        ///
        /// @param region The genomic (0-indexed) region to be cached.
        void cacheSequence( const caller::Region & region );

        const FastaIndex & indexFile() const { return m_indexFile; }

    private:
        void getPaddedSequenceFromFile( std::string * seq, const caller::Region & region ) const;

        /// Compute and return, by reference, the start and end positions, in the FASTA file,
        /// of the specified sequence.
        ///
        /// @param region The genomic region of the sequence
        /// @return The start and end positions in the file.
        std::pair< int64_t, int64_t > computeStartAndEndPosOfSequenceInFile( const caller::Region & region ) const;

    private:
        utils::referenceSequencePtr_t m_cacheSequence;
        mutable ifstream m_file;

        FastaIndex m_indexFile;

        utils::timerPtr_t m_timer;
    };

    using fastaPtr_t = std::shared_ptr< FastaFile >;
}
}

#endif
