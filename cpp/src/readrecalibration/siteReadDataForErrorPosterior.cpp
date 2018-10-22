// All content Copyright (C) 2018 Genomics plc
#include "readrecalibration/siteReadDataForErrorPosterior.hpp"

namespace wecall
{
namespace corrector
{
    void SiteReadDataForErrorPosterior::calculateEmissionProbability(
        const SiteKmerDistribution & siteKmerDistribution,
        ErrorCorrectionParameters const & errorCorrectionParameters )
    {
        const auto & extRefKmer = siteKmerDistribution.getExtReferenceKmer();
        const auto oneOverK = 1.0 / kmerSize;

        double probTrue = 0.0;
        double probError = 0.0;

        const auto & kmerPrior = siteKmerDistribution.kmerPrior();

        for ( const auto & kmerPriorPair : kmerPrior )
        {
            const auto & kmer = kmerPriorPair.first;

            // pad the read kmer with padding from ref kmer at both ends
            // this is done because padding has little effect
            extKmer_t< kmerSize, padding > extTrueKmer;

            for ( std::size_t i = 0; i < padding; ++i )
            {
                extTrueKmer[i] = extRefKmer[i];
                extTrueKmer[i + kmerSize + padding] = extRefKmer[i + kmerSize + padding];
            }

            double probMismatchTrue = 1.0;
            double probMismatchError = 1.0;

            for ( std::size_t i = 0; i < kmerSize; ++i )
            {
                const auto readQualityProb = std::min( 0.75, phred_to_p( m_readQualityKmer[i] ) );

                const auto thisKmerBase_i = kmer[i];
                const auto thisReadBase_i = m_readKmer[i];

                extTrueKmer[i + padding] = thisKmerBase_i;

                if ( thisKmerBase_i == thisReadBase_i )
                {
                    probMismatchTrue *= ( 1.0 - readQualityProb );
                }
                else
                {
                    probMismatchTrue *= ( 1.0 / 3 ) * readQualityProb;
                }
            }

            for ( std::size_t i = 0; i < kmerSize; ++i )
            {
                double probDiff = 0.0;
                const auto thisKmerBase_i = kmer[i];
                const auto thisReadBase_i = m_readKmer[i];

                // with probability pm, a simple match
                if ( thisKmerBase_i == thisReadBase_i )
                {
                    probDiff = errorCorrectionParameters.pm;
                }

                // now deal with the mutation diploid (prob. 1-pm).
                // with probability ps, copy from a distance distribution
                for ( const auto & weightedDistance : errorCorrectionParameters.distances )
                {
                    if ( thisReadBase_i == extTrueKmer[i + padding + weightedDistance.first] )
                    {
                        probDiff += ( 1 - errorCorrectionParameters.pm ) * errorCorrectionParameters.ps *
                                    weightedDistance.second;
                    }
                }

                // with probability of 1-ps, a random mismatch, uniform over 3 possibilities
                if ( thisKmerBase_i != thisReadBase_i )
                {
                    probDiff +=
                        ( 1 - errorCorrectionParameters.pm ) * ( 1 - errorCorrectionParameters.ps ) * ( 1.0 / 3 );
                }

                probMismatchError *= probDiff;
            }

            // true diploid
            probTrue += kmerPriorPair.second * probMismatchTrue;

            // error diploid: include prosisbility of kmer being actually true (and mismatches following Q scores)
            probError += kmerPriorPair.second * ( errorCorrectionParameters.ptrue * probMismatchTrue +
                                                  ( 1 - errorCorrectionParameters.ptrue ) * probMismatchError );

            // WECALL_LOG(SUPER_DEBUG, " Kmer=" << show_string(kmer) << " setPrior=" << kmerPriorPair.second <<
            //                         " probMismTrue=" << probMismatchTrue << " probMismErr=" << probMismatchError);
        }

        // WECALL_LOG(SUPER_DEBUG, " Final emission probs: True=" << probTrue << " Error=" << probError);

        m_emissionProbability.first = fastPow( probTrue, oneOverK );
        m_emissionProbability.second = fastPow( probError, oneOverK );
    }
}
}
