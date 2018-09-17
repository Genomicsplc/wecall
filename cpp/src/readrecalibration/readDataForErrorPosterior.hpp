// All content Copyright (C) 2018 Genomics plc
#ifndef READ_DATA_FOR_ERROR_POSTERIOR_HPP
#define READ_DATA_FOR_ERROR_POSTERIOR_HPP

#include <utility>
#include <vector>

#include "readrecalibration/siteReadDataForErrorPosterior.hpp"
#include "readrecalibration/kmerDistribution.hpp"
#include "readrecalibration/errorCorrectionParameters.hpp"

namespace echidna
{
namespace corrector
{
    class ReadDataForErrorPosterior
    {
    public:
        ReadDataForErrorPosterior( echidna::io::readPtr_t read );

        void calculateProbabilities( const KmerDistribution & kmerDistribution,
                                     const ErrorCorrectionParameters & errorCorrectionParameters );

        void runHmm( KmerDistribution & kmerDistribution );

        /**
         * "Recalibrate" (set to 2) the quality scores where the read is deemed to be in error mode
         */
        void recalibrateRead();

    private:
        double pError( int pos, const KmerDistribution & kmerDistribution ) const
        {
            return kmerDistribution.getSiteKmerDistribution( pos ).pError( m_isForward );
        }

    private:
        echidna::io::readPtr_t m_readPtr;

        std::vector< SiteReadDataForErrorPosterior > m_readData;
        bool m_isForward = false;
        std::size_t m_start = 0;
    };
}  // namespace corrector
}  // namespace echidna

#endif
