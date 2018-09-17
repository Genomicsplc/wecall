// All content Copyright (C) 2018 Genomics plc
#include "utils/logging.hpp"
#include "stats/models.hpp"

#include "stats/functions.hpp"
#include "common.hpp"

namespace echidna
{
namespace stats
{
    double probAllReadsFromSingleCopyOfDiploid( int nReads )
    {
        return betaBinomialCDFForReferenceCalls( nReads, 20.0 );
    }

    double probSufficientVarCoverageToSupportHet( int totalCoverage, int varCoverage )
    {
        if ( totalCoverage == 0 or varCoverage >= ( totalCoverage / 2 ) )
        {
            return std::numeric_limits< double >::quiet_NaN();
        }
        return betaBinomialCDF( varCoverage, totalCoverage, constants::alleleBiasFilterAlpha,
                                constants::alleleBiasFilterBeta );
    }

    double probVarSupportNotBiasedByStrand( int totalFwdCoverage,
                                            int varFwdCoverage,
                                            int totalRevCoverage,
                                            int varRevCoverage )
    {
        const int varCoverage = varFwdCoverage + varRevCoverage;

        if ( varCoverage == 0 or totalRevCoverage == 0 or totalFwdCoverage == 0 )
        {
            return std::numeric_limits< double >::quiet_NaN();
        }
        else
        {
            const double fwdAlpha = constants::strandBiasFilterBeta * totalFwdCoverage / totalRevCoverage;
            const double revAlpha = constants::strandBiasFilterBeta * totalRevCoverage / totalFwdCoverage;

            const auto fwdBias =
                betaBinomialCDF( varFwdCoverage, varCoverage, fwdAlpha, constants::strandBiasFilterBeta );
            const auto revBias =
                betaBinomialCDF( varRevCoverage, varCoverage, revAlpha, constants::strandBiasFilterBeta );

            return std::min( fwdBias, revBias );
        }
    }
}
}
