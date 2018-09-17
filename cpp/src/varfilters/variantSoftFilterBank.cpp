// All content Copyright (C) 2018 Genomics plc
#include "variantSoftFilterBank.hpp"

namespace echidna
{
namespace varfilters
{
    VariantSoftFilterBank::VariantSoftFilterBank( std::vector< std::string > varFilterIDs,
                                                  double alleleBiasThreshP,
                                                  double strandBiasThreshP,
                                                  double allelePlusStrandBiasThreshP,
                                                  phred_t minRootMeanSquareMappingQ,
                                                  double minSNPQOverDepth,
                                                  double minINDELQOverDepth,
                                                  phred_t minBadReadsScore,
                                                  phred_t minCallQual )
    {
        // initialise soft filter vector
        for ( auto filterID : varFilterIDs )
        {
            if ( filterID == vcf::filter::AB_key )
            {
                m_variantFilters.insert( std::make_shared< ABFilter >( alleleBiasThreshP ) );
            }
            else if ( filterID == vcf::filter::SB_key )
            {
                m_variantFilters.insert( std::make_shared< SBFilter >( strandBiasThreshP ) );
            }
            else if ( filterID == vcf::filter::AB_plus_SB_key )
            {
                m_variantFilters.insert( std::make_shared< ABPlusSBFilter >( allelePlusStrandBiasThreshP ) );
            }
            else if ( filterID == vcf::filter::MQ_key )
            {
                m_variantFilters.insert( std::make_shared< MQFilter >( minRootMeanSquareMappingQ ) );
            }
            else if ( filterID == vcf::filter::QD_key )
            {
                m_variantFilters.insert( std::make_shared< QDFilter >( minSNPQOverDepth, minINDELQOverDepth ) );
            }
            else if ( filterID == vcf::filter::BR_key )
            {
                m_variantFilters.insert( std::make_shared< BRFilter >( minBadReadsScore ) );
            }
            else if ( filterID == vcf::filter::LQ_key )
            {
                m_variantFilters.insert( std::make_shared< LQFilter >( minCallQual ) );
            }
            else
            {
                throw utils::echidna_exception( "Invalid variant filter ID: " + filterID );
            }
        }
    }

    std::vector< vcf::FilterDesc > VariantSoftFilterBank::getFilterDescs() const
    {
        std::vector< vcf::FilterDesc > filterDescs;

        for ( const auto & filter : m_variantFilters )
        {
            filterDescs.push_back( filter->getVCFFilterDesc() );
        }

        return filterDescs;
    }

    void VariantSoftFilterBank::applyFilterAnnotation( caller::callVector_t & callSet )
    {
        for ( caller::Call & call : callSet )
        {
            if ( not call.isRefCall() )
            {
                for ( auto filter : m_variantFilters )
                {
                    if ( filter->catches( call ) )
                    {
                        call.filters.insert( filter->getID() );
                    }
                }
            }
        }
    }
}
}
