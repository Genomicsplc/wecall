// All content Copyright (C) 2018 Genomics plc
#ifndef IO_VCF_WRITER_HPP
#define IO_VCF_WRITER_HPP

#include "common.hpp"
#include "vcf/header.hpp"
#include "vcf/record.hpp"
#include "variant/genotype.hpp"
#include "caller/annotation.hpp"
#include "variant/type/variant.hpp"
#include "caller/callSet.hpp"
#include "caller/params.hpp"
#include "io/fastaFile.hpp"
#include "varfilters/filter.hpp"
#include "utils/timer.hpp"

#include <fstream>

namespace echidna
{
namespace io
{
    using variant::genotypePtr_t;
    using variant::varPtr_t;
    using caller::Annotation;

    /// An actor class that serves as an intermediary in the writing out of VCF output.
    class VCFWriter
    {
    public:
        /// Performs file open and sets flags and filters that are consistent through the lifetime of the instance.
        VCFWriter( const std::string & outputFilename, bool outputRefCalls, bool outputPhasedGenotypes );

        /// Destructor to close the open file
        ~VCFWriter();

        /// Creates, and writes out a VCF header.
        ///
        /// @param applicationParams Application parameters as set by user
        /// @param ref The reference file name
        /// @param sampleNames The names of the samples being processed in this job.
        void writeHeader( const std::string & userSpecifiedFormat,
                          const caller::params::Application & applicationParams,
                          const std::string & ref,
                          const std::vector< std::string > & sampleNames,
                          std::vector< vcf::FilterDesc > filterDescs,
                          std::vector< caller::Region > contigs );

        /// Sets state holding the current region being processed - that state being used
        /// by subsequent calls to writeCallSet().
        ///
        /// @param contig Name of the contig (chromosome) being processed.
        void contig( const std::string & contig );

        /// Writes out a set of calls within the current region.
        ///
        /// @param refFile Location of reference genome (for adding indel reference anchor).
        /// @param calls Set of calls to be output.
        void writeCallSet( io::FastaFile & refFile, const caller::callVector_t & calls );

        /// Compiles the REF and ALT strings for one or more calls to be output to a single VCF record.
        ///
        /// @param itB Iterator pointing to the first call to be output.
        /// @param itE Iterator pointing to just past the last call to be output.
        /// @param refFile Location of reference genome (for adding indel reference anchor).
        /// @return A pair containing a single REF string and potentially multiple ALT strings.
        static const std::pair< utils::ReferenceSequence, std::string > compileRefsAndAlts(
            variant::varPtr_t var,
            const io::FastaFile & refFile );

        /// Compiles the annotations of one or more calls for output as INFO on a VCF record.
        ///
        /// @param itB Iterator pointing to the first call to be output.
        /// @param itE Iterator pointing to just past the last call to be output.
        /// @return A structure representing the INFO field of a VCF record.
        static const vcf::Info compileInfo( callIt_t it );

        /// Compiles the sample annotations of one or more calls for output as FORMAT on a VCF record.
        ///
        /// @param itB Iterator pointing to the first call to be output.
        /// @param itE Iterator pointing to just past the last call to be output.
        /// @return A structure representing the FORMAT field of a VCF record.
        static const vcf::SampleInfo compileSampleInfo( callIt_t it, const bool outputPhasedGenotypes );

        /// Compiles the mapping between call and genotype for one or more calls at the same locus.
        ///
        /// @param itB Iterator pointing to the first call to be output.
        /// @param itE Iterator pointing to just past the last call to be output.
        /// @return A structure representing the FORMAT field of a VCF record populated only with the GT field.
        static const std::vector< vcf::SampleInfoValues > compileGenVarMap( callIt_t it,
                                                                            const bool outputPhasedGenotypes );

    private:
        bool m_headerWritten;
        const bool m_outputRefCalls;
        const bool m_outputPhasedGenotypes;
        std::ofstream m_file;

        std::string m_contig;

        utils::timerPtr_t m_timer;
    };
}
}

#endif
