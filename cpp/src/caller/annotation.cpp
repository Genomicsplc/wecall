// All content Copyright (C) 2018 Genomics plc
#include "caller/annotation.hpp"

#include <sstream>
#include <iomanip>
#include "stats/functions.hpp"
#include "vcf/field.hpp"

namespace wecall
{
namespace caller
{
    template <>
    std::vector< std::string > TypedAnnotation< double >::getValues() const
    {
        return {serialise_double( data, def.type )};
    }

    template <>
    std::vector< std::string > TypedAnnotation< std::vector< double > >::getValues() const
    {
        std::vector< std::string > formattedValues;
        for ( auto const & value : data )
        {
            formattedValues.push_back( serialise_double( value, def.type ) );
        }
        return formattedValues;
    }

    template <>
    std::vector< std::string > TypedAnnotation< std::vector< int64_t > >::getValues() const
    {
        std::vector< std::string > formattedValues;
        for ( auto const & value : data )
        {
            formattedValues.push_back( serialise_integral_type( value ) );
        }
        return formattedValues;
    }

    std::string serialise_phred( const phred_t data ) { return serialise_double( data, Annotation::Type::PHRED_VAL ); }

    std::string serialise_double( const double data, const Annotation::Type type )
    {

        if ( std::isnan( data ) or ( type == Annotation::Type::PHRED_VAL and data < 0.0 ) )
        {
            return constants::vcfUnknownValue;
        }
        else if ( type == Annotation::Type::LOG_VAL )
        {
            std::stringstream stringstream;
            stringstream << std::setprecision( 2 ) << std::fixed << data;
            return stringstream.str();
        }
        else if ( type == Annotation::Type::DOUBLE_VAL )
        {
            std::stringstream stringstream;
            stringstream << std::setprecision( 4 ) << data;
            return stringstream.str();
        }
        else if ( type == Annotation::Type::PHRED_VAL )
        {
            return std::to_string( stats::roundPhred( data ) );
        }
        else
        {
            return std::to_string( data );
        }
    }

    std::string serialise_integral_type( const int64_t data ) { return std::to_string( data ); }

    template <>
    std::vector< std::string > TypedAnnotation< int64_t >::getValues() const
    {
        return {serialise_integral_type( data )};
    }

    template <>
    std::vector< std::string > TypedAnnotation< std::vector< std::string > >::getValues() const
    {
        return data;
    }

    template <>
    std::vector< std::string > TypedAnnotation< bool >::getValues() const
    {
        return {std::to_string( data )};
    }

    // Call annotations...
    const Annotation::Def< int64_t > Annotation::BEG =
        Annotation::define< int64_t >( vcf::info::BEG_key, Annotation::COUNT );
    const Annotation::Def< int64_t > Annotation::END =
        Annotation::define< int64_t >( vcf::info::END_key, Annotation::COUNT );
    const Annotation::Def< int64_t > Annotation::LEN =
        Annotation::define< int64_t >( vcf::info::LEN_key, Annotation::COUNT );
    const Annotation::Def< phred_t > Annotation::PP =
        Annotation::define< phred_t >( vcf::info::PP_key, Annotation::PHRED_VAL );
    const Annotation::Def< int64_t > Annotation::DP =
        Annotation::define< int64_t >( vcf::info::DP_key, Annotation::COUNT );
    const Annotation::Def< int64_t > Annotation::DPR =
        Annotation::define< int64_t >( vcf::info::DPR_key, Annotation::COUNT );
    const Annotation::Def< int64_t > Annotation::DPF =
        Annotation::define< int64_t >( vcf::info::DPF_key, Annotation::COUNT );
    const Annotation::Def< int64_t > Annotation::VC =
        Annotation::define< int64_t >( vcf::info::VC_key, Annotation::COUNT );
    const Annotation::Def< int64_t > Annotation::VCR =
        Annotation::define< int64_t >( vcf::info::VCR_key, Annotation::COUNT );
    const Annotation::Def< int64_t > Annotation::VCF =
        Annotation::define< int64_t >( vcf::info::VCF_key, Annotation::COUNT );
    const Annotation::Def< phred_t > Annotation::BR =
        Annotation::define< phred_t >( vcf::info::BR_key, Annotation::PHRED_VAL );
    const Annotation::Def< double > Annotation::ABPV =
        Annotation::define< double >( vcf::info::ABPV_key, Annotation::DOUBLE_VAL );
    const Annotation::Def< double > Annotation::SBPV =
        Annotation::define< double >( vcf::info::SBPV_key, Annotation::DOUBLE_VAL );
    const Annotation::Def< double > Annotation::MQ =
        Annotation::define< double >( vcf::info::MQ_key, Annotation::DOUBLE_VAL );
    const Annotation::Def< double > Annotation::QD =
        Annotation::define< double >( vcf::info::QD_key, Annotation::DOUBLE_VAL );

    // Genotype call annotations...
    const Annotation::Def< std::vector< phred_t > > Annotation::PL =
        Annotation::define< std::vector< phred_t > >( vcf::format::PL_key, Annotation::PHRED_VAL );

    const Annotation::Def< std::vector< int64_t > > Annotation::AD =
        Annotation::define< std::vector< int64_t > >( vcf::format::AD_key, Annotation::COUNT );

    const Annotation::Def< std::vector< double > > Annotation::VAF =
        Annotation::define< std::vector< double > >( vcf::format::VAF_key, Annotation::DOUBLE_VAL );

    const Annotation::Def< phred_t > Annotation::GQ =
        Annotation::define< phred_t >( vcf::format::GQ_key, Annotation::PHRED_VAL );
    const Annotation::Def< int64_t > Annotation::FORMAT_DP =
        Annotation::define< int64_t >( vcf::format::DP_key, Annotation::DOUBLE_VAL );
    const Annotation::Def< int64_t > Annotation::MIN_DP =
        Annotation::define< int64_t >( vcf::format::MIN_DP_key, Annotation::COUNT );
    const Annotation::Def< int64_t > Annotation::PS =
        Annotation::define< int64_t >( vcf::format::PS_key, Annotation::COUNT );
    const Annotation::Def< phred_t > Annotation::PQ =
        Annotation::define< phred_t >( vcf::format::PQ_key, Annotation::PHRED_VAL );
}
}
