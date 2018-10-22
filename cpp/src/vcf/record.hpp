// All content Copyright (C) 2018 Genomics plc
#ifndef VCF_RECORD_HPP
#define VCF_RECORD_HPP

#include "common.hpp"
#include "caller/annotation.hpp"
#include "variant/type/variant.hpp"

#include <ostream>
#include <set>

namespace wecall
{
namespace vcf
{
    using Info = std::vector< std::pair< std::string, std::vector< std::string > > >;
    using SampleInfoFormat = std::vector< std::string >;
    using SampleInfoValues = std::vector< std::vector< std::string > >;
    using SampleInfo = std::pair< SampleInfoFormat, std::vector< SampleInfoValues > >;

    /// An in-memory representation of a VCF record.
    class Record
    {
    public:
        /// Basic constructor from constituent data.
        ///
        /// @param contig Contiguous DNA segment name - typically chromosome name.
        /// @param pos Base position (locus) within the contig (1-indexed).
        /// @param id Variant identifier ("." if unknown).
        /// @param ref Reference sequence at the locus.
        /// @param alts List of alternate allele sequences.
        /// @param qual Quality of the call expressed as a Phred score.
        /// @param filters List of filters that this call fails, or "PASS" if call is good.
        /// @param info Set of call annotations.
        /// @param sampleInfo Sets of sample-specific call annotations.
        Record( std::string contig,
                std::size_t pos,
                std::set< std::string > ids,
                std::string ref,
                std::vector< std::string > alts,
                phred_t qual,
                std::set< std::string > filters,
                Info info,
                SampleInfo sampleInfo );

        std::vector< variant::varPtr_t > getVariants( const utils::referenceSequencePtr_t & referenceSequence ) const;

        /// Writes the VCF record to the output stream.
        ///
        /// @param out Output stream
        /// @param record Record to be output
        /// @return Output stream
        friend std::ostream & operator<<( std::ostream & out, const Record & record );

        const std::string m_contig;
        const std::size_t m_pos;
        const std::set< std::string > m_ids;
        const std::string m_ref;
        const std::vector< std::string > m_alts;
        const double m_qual;
        const std::set< std::string > m_filters;
        const Info m_info;
        const SampleInfo m_sampleInfo;
    };
}
}

#endif
