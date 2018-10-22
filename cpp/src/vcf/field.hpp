// All content Copyright (C) 2018 Genomics plc
#ifndef VCF_FIELD_HPP
#define VCF_FIELD_HPP

#include <string>
#include <vector>
#include <iostream>

#include "common.hpp"

namespace wecall
{
namespace vcf
{

    namespace info
    {
        const static std::string VC_key = "VC";
        const static std::string VCF_key = "VCF";
        const static std::string VCR_key = "VCR";

        const static std::string DP_key = "DP";
        const static std::string DPR_key = "DPR";
        const static std::string DPF_key = "DPF";

        const static std::string PP_key = "PP";

        const static std::string ABPV_key = "ABPV";
        const static std::string SBPV_key = "SBPV";

        const static std::string MQ_key = "MQ";
        const static std::string QD_key = "QD";
        const static std::string BR_key = "BR";

        const static std::string BEG_key = "BEG";
        const static std::string END_key = "END";
        const static std::string LEN_key = "LEN";

        std::vector< std::string > getVCFKeys( bool outputRefCalls );
    }

    namespace filter
    {
        const static std::string AB_key = "AB";
        const static std::string SB_key = "SB";
        const static std::string AB_plus_SB_key = "AB+SB";
        const static std::string BR_key = "BR";
        const static std::string MQ_key = "MQ";
        const static std::string QD_key = "QD";
        const static std::string NC_key = "NC";
        const static std::string LQ_key = "LQ";

        const static std::vector< std::string > VCFKeys =
            {AB_key, SB_key, AB_plus_SB_key, MQ_key, QD_key, BR_key, LQ_key};
    }

    namespace format
    {
        const static std::string GT_key = "GT";
        const static std::string AD_key = "AD";
        const static std::string DP_key = "DP";
        const static std::string VAF_key = "VAF";
        const static std::string MIN_DP_key = "MIN_DP";
        const static std::string GQ_key = "GQ";
        const static std::string PL_key = "PL";
        const static std::string PS_key = "PS";
        const static std::string PQ_key = "PQ";

        std::vector< std::string > getVCFKeys( bool outputPhasedGenotypes, bool outputRefCalls );
    }

    /// An in-memory representation of a VCF field (annotation) definition.
    class Field
    {
    public:
        /// Look-up the definition of a VCF INFO field given its ID.
        ///
        /// @param ID The field ID.
        /// @param source The source name of the application generating the vcf output
        /// @param version The version number of the application generating the vcf output
        /// @return The field definition.
        static Field infoFieldFromID( const std::string & ID,
                                      const std::string & source = "",
                                      const std::string & version = "" );

        /// Look-up the definition of a VCF FORMAT field given its ID.
        ///
        /// @param ID The field ID.
        /// @param source The source name of the application generating the vcf output
        /// @param version The version number of the application generating the vcf output
        /// @return The field definition.
        static Field formatFieldFromID( const std::string & ID,
                                        const std::string & source = "",
                                        const std::string & version = "" );

        /// Options for the cardinality of the field.
        enum Cardinality
        {
            UNKNOWN,
            FIXED_NO,
            ALLELE,
            ALT_ALLELE,
            GENOTYPE
        };

        /// Options for the field type.
        enum Type
        {
            INTEGER,
            FLOAT,
            FLAG,
            CHARACTER,
            STRING
        };

        /// Basic constructor from constituent data
        ///
        /// @param id Field ID.
        /// @param cardinality Together with number, maps to VCF INFO/FORMAT Number attribute.
        /// @param number Together with cardinality, maps to VCF INFO/FORMAT Number attribute.
        /// @param type Field data type.
        /// @param description Field description.
        /// @param source Field source (typically EchiDNA). Only written onto INFO output
        /// @param version Field version. Only written onto INFO output
        Field( std::string id,
               Cardinality cardinality,
               std::size_t number,
               Type type,
               std::string description,
               std::string source,
               std::string version );

        /// Writes the field to the output stream in the form of a VCF INFO/FORMAT header line.
        ///
        /// @param out Output stream
        /// @param field Field to be output
        /// @return Output stream
        friend std::ostream & operator<<( std::ostream & out, const Field & field );

        /// Tests whether or not a field is single-valued.
        ///
        /// @return True if field is single-valued.
        bool isSingleValued() const { return ( m_cardinality == FIXED_NO and m_number == 1 ); }

        /// Tests whether or not a field has one value per ALT allele.
        ///
        /// @return True if field has one value per ALT allele
        bool hasAltCardinality() const { return m_cardinality == ALT_ALLELE; }

        const std::string m_id;
        const Cardinality m_cardinality;
        const std::size_t m_number;
        const Type m_type;
        const std::string m_description;
        const std::string m_source;
        const std::string m_version;
    };
}
}

#endif
