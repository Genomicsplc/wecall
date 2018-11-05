// All content Copyright (C) 2018 Genomics plc
#ifndef VARFILTERS_VARIANTFILTERBANK_HPP
#define VARFILTERS_VARIANTFILTERBANK_HPP

#include "varfilters/filter.hpp"
#include "caller/region.hpp"
#include "caller/callSet.hpp"

namespace wecall
{
namespace varfilters
{
    using varFilterId_t = std::string;

    class VariantSoftFilterBank
    {
    public:
        VariantSoftFilterBank( std::vector< std::string > varFilterIDs,
                               double alleleBiasThreshP,
                               double strandBiasThreshP,
                               double allelePlusStrandBiasThreshP,
                               phred_t minRootMeanSquareMappingQ,
                               double minSNPQOverDepth,
                               double minINDELQOverDepth,
                               phred_t minBadReadsScore,
                               phred_t minCallQual );

        std::vector< vcf::FilterDesc > getFilterDescs() const;

        void applyFilterAnnotation( caller::callVector_t & callSet );

    private:
        std::set< varfilters::FilterPtr_t, FilterPtrComp > m_variantFilters;
        std::vector< vcf::FilterDesc > m_systematicFilterDescs;
    };
}
}

#endif  // VARFILTERS_VARIANTFILTERBANK_HPP
