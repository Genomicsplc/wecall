// All content Copyright (C) 2018 Genomics plc
#ifndef VARIANT_NORMALIZER_HPP
#define VARIANT_NORMALIZER_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/optional/optional.hpp>

#include "utils/referenceSequence.hpp"
#include "utils/NeedlemanWunsch.hpp"
#include "variant/type/variant.hpp"

namespace wecall
{
namespace variant
{
    class VariantNormalizer
    {
        using scoreMatrix_t = boost::numeric::ublas::matrix< int32_t >;

        using bitType_t = int64_t;
        using traceMatrix_t = boost::numeric::ublas::matrix< bitType_t >;

    public:
        VariantNormalizer( const utils::referenceSequencePtr_t & referenceSequence )
            : m_referenceSequence( referenceSequence )
        {
        }

        variantSet_t getNormalized( const caller::Region & region,
                                    const utils::BasePairSequence & sequence,
                                    const boost::optional< variantSet_t > unnormalizedVariants ) const;

    private:
        const utils::referenceSequencePtr_t m_referenceSequence;
        const uint64_t m_maxMatrixSize = 50000;
    };
}
}
#endif
