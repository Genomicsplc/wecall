// All content Copyright (C) 2018 Genomics plc
#ifndef BAMFILE_HPP
#define BAMFILE_HPP

#include <string>
#include <ctype.h>
#include <assert.h>

// Declare the samtools headers as extern
extern "C" {
#include "samtools/sam.h"
#include "samtools/bam.h"
#include "samtools/khash.h"
#include "samtools/ksort.h"
#include "tabix/bam_endian.h"
#include "samtools/knetfile.h"
}

#include <map>
#include <vector>
#include <samtools/sam.h>

#include "utils/logging.hpp"
#include "utils/timer.hpp"
#include "caller/region.hpp"
#include "io/readDataSet.hpp"
#include "io/read.hpp"
#include "bamFileIterator.hpp"
#include "io/pysam.hpp"

namespace echidna
{
namespace io
{
    using readVector_t = std::vector< echidna::io::readPtr_t >;
    using readMap_t = std::map< std::string, readVector_t >;

    /// Concrete data source class to represent a BAM file. BAM is a compressed (using bgzip)
    /// binary format for storing read data, typically (but not neccessarily) aligned to a
    /// reference genome. Here we use Samtools (http://samtools.sourceforge.net/) to provde
    /// random access to the read data.
    class BamFile
    {
    public:
        /// Construct a BamFile instance. If the file exists, it wil be opened for
        /// reading.
        ///
        /// @param fileName The name of the BAM file in the filesystem. Must be readable
        explicit BamFile( std::string fileName );

        /// Destructor. Free up memory and close file
        ~BamFile();

        bool isOpen() const { return m_samFile != 0; }

        std::string getContigName( const std::size_t tid ) const;

        std::vector< std::string > getSampleNames() const;

        /// Get the sample names mapped to read-group IDs
        ///
        /// @return A map of read-group IDs to sample names
        std::map< std::string, std::string > getSamplesByID() const;

        bamFileIteratorPtr_t readRegion( const caller::Region & blockRegion,
                                         utils::referenceSequencePtr_t m_refSequence );

        bam_header_t bamHeader() const { return *( m_samFile->header ); }

    private:
        std::string m_fileName;  /// The name of the file
        samfile_t * m_samFile;   /// A pointer to the Samtools data-structure
        bam_index_t * m_index;   /// Pointer to the Samtools index structure
        utils::timerPtr_t m_timer;
        std::map< std::string, std::string > m_samplesByID;

        std::string default_sample_name() const;
    };

    using bamFilePtr_t = std::shared_ptr< BamFile >;
}
}

#endif
