// All content Copyright (C) 2018 Genomics plc
#include "vcf/field.hpp"
#include "utils/exceptions.hpp"

#include <stdexcept>

namespace wecall
{
namespace vcf
{
    namespace info
    {
        std::vector< std::string > getVCFKeys( bool outputRefCalls )
        {
            std::vector< std::string > vcfKeys = {
                PP_key, DP_key, DPR_key, DPF_key, VC_key, VCR_key, VCF_key, MQ_key, ABPV_key, SBPV_key, QD_key, BR_key};
            if ( outputRefCalls )
            {
                const std::vector< std::string > refCallOnlyKeys = {BEG_key, END_key, LEN_key};
                vcfKeys.insert( vcfKeys.end(), refCallOnlyKeys.cbegin(), refCallOnlyKeys.cend() );
            }
            return vcfKeys;
        }
    }

    namespace format
    {
        std::vector< std::string > getVCFKeys( bool outputPhasedGenotypes, bool outputRefCalls )
        {
            std::vector< std::string > vcfKeys = {GT_key, AD_key, DP_key, GQ_key, PL_key, VAF_key};
            if ( outputPhasedGenotypes )
            {
                const std::vector< std::string > phasedOnlyKeys = {PS_key, PQ_key};
                vcfKeys.insert( vcfKeys.end(), phasedOnlyKeys.cbegin(), phasedOnlyKeys.cend() );
            }
            if ( outputRefCalls )
            {
                vcfKeys.push_back( MIN_DP_key );
            }
            return vcfKeys;
        }
    }

    // Keys reused in multiple descriptions.

    // INFO fields...
    // Include in INFO fields the application name and version (e.g. EchiDNA, 1.0)
    Field Field::infoFieldFromID( const std::string & ID, const std::string & source, const std::string & version )
    {
        if ( ID == info::PP_key )
        {
            return Field( ID, Field::ALT_ALLELE, 1, Field::INTEGER,
                          "Posterior probability (phred scaled) that this variant does not segregate.", source,
                          version );
        }
        else if ( ID == info::DP_key )
        {
            return Field( ID, Field::FIXED_NO, 1, Field::INTEGER, "Total depth of read coverage at this locus.", source,
                          version );
        }
        else if ( ID == info::DPF_key )
        {
            return Field( ID, Field::FIXED_NO, 1, Field::INTEGER,
                          "Total probabilistic depth of forward read coverage at this locus (sum of probabilities of "
                          "each read supporting the variant).",
                          source, version );
        }
        else if ( ID == info::DPR_key )
        {
            return Field( ID, Field::FIXED_NO, 1, Field::INTEGER,
                          "Total probabilistic depth of reverse read coverage at this locus (sum of probabilities of "
                          "each read supporting the variant).",
                          source, version );
        }
        else if ( ID == info::VC_key )
        {
            return Field( ID, Field::ALT_ALLELE, 1, Field::INTEGER,
                          "Total probabilistic number of reads supporting each alternative allele (sum of "
                          "probabilities of each read supporting the allele).",
                          source, version );
        }
        else if ( ID == info::VCF_key )
        {
            return Field( ID, Field::ALT_ALLELE, 1, Field::INTEGER,
                          "Total probabilistic number of forward reads supporting each alternative allele (sum of "
                          "probabilities of each read supporting the allele).",
                          source, version );
        }
        else if ( ID == info::VCR_key )
        {
            return Field( ID, Field::ALT_ALLELE, 1, Field::INTEGER,
                          "Total probabilistic number of reverse reads supporting each alternative allele (sum of "
                          "probabilities of each read supporting the allele).",
                          source, version );
        }
        else if ( ID == info::ABPV_key )
        {
            return Field( ID, Field::ALT_ALLELE, 1, Field::FLOAT,
                          "Allele bias P-value; probability that fraction of reads supporting alt allele (" +
                              info::VC_key + ") amongst read depth (" + info::DP_key +
                              ") is more extreme than expected assuming a beta-binomial distribution.",
                          source, version );
        }
        else if ( ID == info::SBPV_key )
        {
            return Field( ID, Field::ALT_ALLELE, 1, Field::FLOAT,
                          "Strand bias P-value; probability that the fraction of forward reads (" + info::VCF_key +
                              ") amongst reads supporting alt allele (" + info::VC_key +
                              ") is more extreme than expected assuming a beta-binomial distribution.",
                          source, version );
        }
        else if ( ID == info::BEG_key )
        {
            return Field( ID, Field::FIXED_NO, 1, Field::INTEGER, "Start position of reference call block.", source,
                          version );
        }
        else if ( ID == info::END_key )
        {
            return Field( ID, Field::FIXED_NO, 1, Field::INTEGER, "End position of reference call block (inclusive).",
                          source, version );
        }
        else if ( ID == info::LEN_key )
        {
            return Field( ID, Field::FIXED_NO, 1, Field::INTEGER, "Length of reference call block.", source, version );
        }
        else if ( ID == info::MQ_key )
        {
            return Field( ID, Field::ALT_ALLELE, 1, Field::FLOAT,
                          "Root mean square of mapping quality of reads supporting each alternative allele.", source,
                          version );
        }
        else if ( ID == info::BR_key )
        {
            return Field( ID, Field::ALT_ALLELE, 1, Field::FLOAT,
                          "The median of the per-read min base quality (within a interval of the locus) taken over "
                          "reads supporting each allele.",
                          source, version );
        }
        else if ( ID == info::QD_key )
        {
            return Field( ID, Field::ALT_ALLELE, 1, Field::FLOAT,
                          "Ratio of phred-scaled posterior probability (" + info::PP_key +
                              ") to number of supporting reads for each allele (" + info::VC_key + ").",
                          source, version );
        }
        else
        {
            throw utils::wecall_exception( "Invalid field ID: " + ID );
        }
    }

    // FORMAT fields...
    // Include in FORMAT fields the application name and version for FORMAT fields
    Field Field::formatFieldFromID( const std::string & ID, const std::string & source, const std::string & version )
    {
        if ( ID == format::GT_key )
        {
            return Field( ID, Field::FIXED_NO, 1, Field::STRING,
                          "Genotypes of reference and alternative alleles in order listed.", source, version );
        }
        else if ( ID == format::PL_key )
        {
            // TODO(ES): This number 3 is a magic number! It is ploidy + 1.
            return Field( ID, Field::GENOTYPE, 3, Field::INTEGER,
                          "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification.",
                          source, version );
        }
        else if ( ID == format::AD_key )
        {
            return Field( ID, Field::UNKNOWN, 3, Field::INTEGER,
                          "Probabilistic allelic depths for the ref and alt alleles in the order listed (i.e. INFO::" +
                              info::VC_key + " split out by sample).",
                          source, version );
        }
        else if ( ID == format::GQ_key )
        {
            return Field(
                ID, Field::FIXED_NO, 1, Field::INTEGER,
                "Phred-scaled genotype quality (i.e. posterior probability that the genotype call is incorrect).",
                source, version );
        }
        else if ( ID == format::DP_key )
        {
            return Field( ID, Field::FIXED_NO, 1, Field::INTEGER,
                          "Number of reads overlapping the variant site (i.e. INFO::" + info::DP_key +
                              " split out by sample). For reference calls the average depth (rounded to the nearest "
                              "integer) over the region is reported.",
                          source, version );
        }
        else if ( ID == format::MIN_DP_key )
        {
            return Field( ID, Field::FIXED_NO, 1, Field::INTEGER,
                          "Minimum read coverage observed within the reference block.", source, version );
        }
        else if ( ID == format::PS_key )
        {
            return Field( ID, Field::FIXED_NO, 1, Field::STRING, "Phase set id.", source, version );
        }
        else if ( ID == format::PQ_key )
        {
            return Field( ID, Field::FIXED_NO, 1, Field::INTEGER,
                          "Phred-scaled phase quality (i.e. posterior probability that the phasing is incorrect).",
                          source, version );
        }
        else if ( ID == format::VAF_key )
        {
            return Field( ID, Field::ALT_ALLELE, 1, Field::FLOAT,
                          "Probabilistic variant allelic frequencies for each alt allele (FORMAT::" + format::AD_key +
                              " / FORMAT::" + format::DP_key + ").",
                          source, version );
        }
        else
        {
            throw utils::wecall_exception( "Invalid field ID: " + ID );
        }
    }

    Field::Field( std::string id,
                  Cardinality cardinality,
                  std::size_t number,
                  Type type,
                  std::string description,
                  std::string source,
                  std::string version )
        : m_id( id ),
          m_cardinality( cardinality ),
          m_number( number ),
          m_type( type ),
          m_description( description ),
          m_source( source ),
          m_version( version )
    {
        // Nothing to do here
    }

    std::ostream & operator<<( std::ostream & out, const Field & field )
    {
        out << "<ID=" << field.m_id;

        out << ",Number=";
        switch ( field.m_cardinality )
        {
        case Field::UNKNOWN:
            out << constants::vcfUnknownValue;
            break;
        case Field::FIXED_NO:
            out << field.m_number;
            break;
        case Field::ALLELE:
            out << "R";
            break;
        case Field::ALT_ALLELE:
            out << "A";
            break;
        case Field::GENOTYPE:
            out << "G";
            break;
        }

        out << ",Type=";
        switch ( field.m_type )
        {
        case Field::INTEGER:
            out << "Integer";
            break;
        case Field::FLOAT:
            out << "Float";
            break;
        case Field::FLAG:
            out << "Flag";
            break;
        case Field::CHARACTER:
            out << "Character";
            break;
        case Field::STRING:
            out << "String";
            break;
        }

        out << ",Description=\"" << field.m_description << "\"";

        // To include application source name and version for INFO fields. Source name and version need to be in
        // quotation marks.
        if ( not field.m_source.empty() and not field.m_version.empty() )
        {
            out << ",Source=\"" << field.m_source << "\"";
            out << ",Version=\"" << field.m_version << "\"";
        }

        out << ">";

        return out;
    }
}
}
