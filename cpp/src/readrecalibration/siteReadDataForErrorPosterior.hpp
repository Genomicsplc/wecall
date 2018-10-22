// All content Copyright (C) 2018 Genomics plc
#ifndef SITE_READ_DATA_FOR_ERROR_POSTERIOR_HPP
#define SITE_READ_DATA_FOR_ERROR_POSTERIOR_HPP

#include <utility>
#include <vector>

#include "readrecalibration/siteKmerDistribution.hpp"
#include "utils/sequence.hpp"

namespace wecall
{
namespace corrector
{

    //
    // This data structure is associated with a particular nucleotide in a read
    //
    class SiteReadDataForErrorPosterior
    {

    public:
        SiteReadDataForErrorPosterior( const wecall::utils::BasePairSequence & sequenceKmer,
                                       const wecall::utils::BasePairSequence & qualityKmer,
                                       std::size_t indexIntoRead,
                                       std::size_t refPos )
        {
            std::copy( sequenceKmer.cbegin(), sequenceKmer.cend(), m_readKmer.data() );
            std::copy( qualityKmer.cbegin(), qualityKmer.cend(), m_readQualityKmer.data() );

            m_indexIntoRead = indexIntoRead;
            m_referencePosition = refPos;
        }

        void calculateEmissionProbability( const SiteKmerDistribution & siteKmerDistribution,
                                           const ErrorCorrectionParameters & errorCorrectionParameters );

        void calculateTransitionProbability( ErrorCorrectionParameters const & errorCorrectionParameters )
        {
            // TODO(sem): calculate error posterior at iteration #2
            m_pErrorForward = errorCorrectionParameters.perr;
            m_pErrorBackward = errorCorrectionParameters.perr;
        }

        std::size_t indexIntoRead() const { return m_indexIntoRead; }
        std::size_t referencePosition() const { return m_referencePosition; }

        double errorPosterior() { return m_errorPosterior; }
        void errorPosterior( double errorPosterior ) { m_errorPosterior = errorPosterior; }

        double errorTransitionProbability() { return m_errorTransitionProbability; }
        void errorTransitionProbability( double errorTransistionProbability )
        {
            m_errorTransitionProbability = errorTransistionProbability;
        }

    public:
        using probabilityVectorOfLengthTwo_t = std::pair< double, double >;
        using partialLikelihood_t = probabilityVectorOfLengthTwo_t;  ///< The likelihood for normal and error states
        using emissionProbability_t = probabilityVectorOfLengthTwo_t;

        emissionProbability_t emissionProbability() const { return m_emissionProbability; }
        partialLikelihood_t & forwardLikelihood() { return m_forwardLikelihood; }

    private:
        // inputs
        double m_pErrorForward = 0;
        double m_pErrorBackward = 0;
        std::size_t m_indexIntoRead = 0;
        std::size_t m_referencePosition = 0;
        emissionProbability_t m_emissionProbability = {0, 0};

        // intermediate results
        partialLikelihood_t m_forwardLikelihood = {0, 0};

        // outputs
        double m_errorPosterior = 0;              ///< Posterior probability of being in the error state
        double m_errorTransitionProbability = 0;  ///< Posterior probability of transition into the error state

        kmer_t< kmerSize > m_readKmer;         ///< Kmer at this read position
        kmer_t< kmerSize > m_readQualityKmer;  ///< Kmer of quality values at this position
    };
}  // namespace corrector
}  // namespace wecall

#endif
