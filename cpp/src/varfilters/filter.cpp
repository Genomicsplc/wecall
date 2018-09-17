// All content Copyright (C) 2018 Genomics plc
#include "varfilters/filter.hpp"

#include <stdexcept>
#include "vcf/filterDescription.hpp"
#include "utils/write.hpp"
#include "vcf/field.hpp"

namespace echidna
{
namespace varfilters
{

    bool filterImpl( double value, double threshold )
    {
        if ( std::isnan( value ) )
        {
            return false;
        }
        else
        {
            return value < threshold;
        }
    }

    Filter::Filter( const std::string & id, const std::string & description ) : m_id( id ), m_description( description )
    {
    }

    vcf::FilterDesc Filter::getVCFFilterDesc() { return vcf::FilterDesc( m_id, m_description ); }

    bool Filter::catches( const caller::Call & call ) { return false; }

    std::ostream & operator<<( std::ostream & out, const Filter & filter )
    {
        out << "<ID=" << filter.m_id << ",Description=\"" << filter.m_description << "\">";

        return out;
    }

    ABFilter::ABFilter( double thresholdP )
        : Filter( vcf::filter::AB_key,
                  "Allele Bias: Indicates lower number of reads supporting variant than expected (any of INFO::" +
                      vcf::info::ABPV_key +
                      " < " +
                      utils::toString( thresholdP ) +
                      ")." ),
          m_thresholdP( thresholdP )
    {
    }

    bool ABFilter::catches( const caller::Call & call )
    {
        return filterImpl( call.getAnnotation( Annotation::ABPV ), m_thresholdP );
    }

    SBFilter::SBFilter( double thresholdP )
        : Filter( vcf::filter::SB_key,
                  "Strand Bias: Indicates imbalance between number of forward and reverse reads supporting variant "
                  "(any of INFO::" +
                      vcf::info::SBPV_key +
                      " < " +
                      utils::toString( thresholdP ) +
                      ")." ),
          m_thresholdP( thresholdP )
    {
    }

    bool SBFilter::catches( const caller::Call & call )
    {
        return filterImpl( call.getAnnotation( Annotation::SBPV ), m_thresholdP );
    }

    ABPlusSBFilter::ABPlusSBFilter( double thresholdP )
        : Filter( vcf::filter::AB_plus_SB_key,
                  "Allele + Strand Bias: Indicates that both the AB and SB filters are close to being triggered "
                  "(any of INFO::" +
                      vcf::info::ABPV_key +
                      " + INFO::" +
                      vcf::info::SBPV_key +
                      " < " +
                      utils::toString( thresholdP ) +
                      ")." ),
          m_thresholdP( thresholdP )
    {
    }

    bool ABPlusSBFilter::catches( const caller::Call & call )
    {
        const double abpv = call.getAnnotation( Annotation::ABPV );
        const double sbpv = call.getAnnotation( Annotation::SBPV );
        return filterImpl( abpv + sbpv, m_thresholdP );
    }

    MQFilter::MQFilter( double thresholdP )
        : Filter( vcf::filter::MQ_key,
                  "low Mapping Quality: Indicates presence of low mapping quality (any of INFO::" + vcf::info::MQ_key +
                      " < " +
                      utils::toString( thresholdP ) +
                      ")." ),
          m_thresholdP( thresholdP )
    {
    }

    bool MQFilter::catches( const caller::Call & call )
    {
        return filterImpl( call.getAnnotation( Annotation::MQ ), m_thresholdP );
    }

    QDFilter::QDFilter( double SNPThresholdP, double INDELThresholdP )
        : Filter( vcf::filter::QD_key,
                  "Quality over Depth: Indicates low quality relative to number of supporting reads (any of INFO::" +
                      vcf::info::QD_key +
                      " < " +
                      utils::toString( INDELThresholdP ) +
                      " for Indels or INFO::" +
                      vcf::info::QD_key +
                      " < " +
                      utils::toString( SNPThresholdP ) +
                      " otherwise)." ),
          m_SNPThresholdP( SNPThresholdP ),
          m_INDELThresholdP( INDELThresholdP )
    {
    }

    bool QDFilter::catches( const caller::Call & call )
    {
        if ( call.var->isIndel() )
        {
            return filterImpl( call.getAnnotation( Annotation::QD ), m_INDELThresholdP );
        }
        else
        {
            return filterImpl( call.getAnnotation( Annotation::QD ), m_SNPThresholdP );
        }
    }

    BRFilter::BRFilter( phred_t thresholdPhred )
        : Filter(
              vcf::filter::BR_key,
              "Bad Reads: Indicates low quality base pairs on reads in the vicinity of variant locus (any of INFO::" +
                  vcf::info::BR_key +
                  " < " +
                  utils::toString( thresholdPhred ) +
                  ")." ),
          m_thresholdPhred( thresholdPhred )
    {
    }

    bool BRFilter::catches( const caller::Call & call )
    {
        return call.getAnnotation( Annotation::BR ) < m_thresholdPhred;
    }

    LQFilter::LQFilter( phred_t thresholdPhred )
        : Filter( vcf::filter::LQ_key,
                  "Low Quality: Indicates a low variant quality (any of INFO::" + vcf::info::PP_key + " < " +
                      utils::toString( thresholdPhred ) +
                      ")." ),
          m_thresholdPhred( thresholdPhred )
    {
    }

    bool LQFilter::catches( const caller::Call & call )
    {
        return call.getAnnotation( Annotation::PP ) < m_thresholdPhred;
    }

    bool FilterPtrComp::operator()( const FilterPtr_t & x, const FilterPtr_t & y ) const
    {
        return x->getID() < y->getID();
    }
}
}
