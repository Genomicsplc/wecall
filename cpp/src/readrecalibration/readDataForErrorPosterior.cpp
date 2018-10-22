// All content Copyright (C) 2018 Genomics plc
#include "readrecalibration/readDataForErrorPosterior.hpp"
#include "alignment/cigarItems.hpp"

namespace wecall
{
namespace corrector
{
    ReadDataForErrorPosterior::ReadDataForErrorPosterior( wecall::io::readPtr_t readPtr ) : m_readPtr( readPtr )
    {
        m_isForward = not m_readPtr->isReverse();

        auto referencePositions = m_readPtr->getReferencePositions();

        for ( size_t indexIntoRead = 0; indexIntoRead < referencePositions.size() - kmerSize + 1; ++indexIntoRead )
        {
            const int refPos = referencePositions[indexIntoRead];

            if ( refPos == alignment::emptyPos )
            {
                continue;
            }
            else
            {
                auto kmer = m_readPtr->sequence().substr( indexIntoRead, kmerSize );
                auto readQualityKmer = m_readPtr->getQualities().substr( indexIntoRead, kmerSize );

                m_readData.emplace_back( kmer, readQualityKmer, indexIntoRead, refPos );
            }
        }
    }

    void ReadDataForErrorPosterior::calculateProbabilities(
        const KmerDistribution & kmerDistribution,
        const ErrorCorrectionParameters & errorCorrectionParameters )
    {
        for ( auto & siteReadData : m_readData )
        {
            int refPos = siteReadData.referencePosition();

            if ( kmerDistribution.start() <= refPos && refPos < kmerDistribution.end() )
            {
                // ECHIDNA_LOG(SUPER_DEBUG,
                //            "Emission distribution for read=" << m_readPtr->getQName() << " start position=" <<
                //            m_readPtr->getStartPos());
                siteReadData.calculateEmissionProbability( kmerDistribution.getSiteKmerDistribution( refPos ),
                                                           errorCorrectionParameters );

                // ECHIDNA_LOG(SUPER_DEBUG,
                //            "Emission transition for read=" << m_readPtr->getQName() << " start position=" <<
                //            m_readPtr->getStartPos());
                siteReadData.calculateTransitionProbability( errorCorrectionParameters );
            }
        }

        m_start = kmerDistribution.start();
    }

    void ReadDataForErrorPosterior::runHmm( KmerDistribution & kmerDistribution )
    {
        if ( m_readData.size() == 0 )
        {
            return;
        }

        int const startIndex = m_isForward ? 0 : m_readData.size() - 1;  // inclusive index
        int const endIndex = m_isForward ? m_readData.size() - 1 : 0;    // inclusive index
        int const direction = m_isForward ? 1 : -1;

        double likelihoodTrue = 1.0;
        double likelihoodError = 0.0;

        for ( int i = startIndex; /* empty condition */; i += direction )
        {
            const auto pErrorForPos = pError( m_readData[i].referencePosition(), kmerDistribution );

            likelihoodError =
                ( likelihoodError + likelihoodTrue * pErrorForPos ) * m_readData[i].emissionProbability().second;
            likelihoodTrue = likelihoodTrue * ( 1 - pErrorForPos ) * m_readData[i].emissionProbability().first;

            m_readData[i].forwardLikelihood() = {likelihoodTrue, likelihoodError};

            if ( i == endIndex )
            {
                break;
            }
        }

        double const likelihoodSum = likelihoodTrue + likelihoodError;
        std::pair< double, double > backwardLikelihood = {1.0 / likelihoodSum, 1.0 / likelihoodSum};

        for ( int i = endIndex; /* empty condition */; i -= direction )
        {
            const auto pErrorForPos = pError( m_readData[i].referencePosition(), kmerDistribution );

            m_readData[i].errorPosterior( backwardLikelihood.second * m_readData[i].forwardLikelihood().second );

            backwardLikelihood.first =
                m_readData[i].emissionProbability().first * ( 1 - pErrorForPos ) * backwardLikelihood.first +
                m_readData[i].emissionProbability().second * pErrorForPos * backwardLikelihood.second;

            backwardLikelihood.second = m_readData[i].emissionProbability().second * backwardLikelihood.second;

            if ( i == startIndex )
            {
                break;
            }
        }

        // calculate the probabilities of moving into the error state
        double currentErrorPosterior = 0.0;
        for ( int i = startIndex; /* empty condition */; i += direction )
        {
            int const pos = m_readData[i].referencePosition();

            double errorProbability = m_readData[i].errorPosterior() - currentErrorPosterior;

            m_readData[i].errorTransitionProbability( errorProbability );

            kmerDistribution.accumulateErrorProbability( pos, m_readData[i].errorTransitionProbability(), m_isForward );

            currentErrorPosterior = m_readData[i].errorPosterior();

            if ( i == endIndex )
            {
                break;
            }
        }
    }

    void ReadDataForErrorPosterior::recalibrateRead()
    {
        // ECHIDNA_LOG(SUPER_DEBUG, "Before recalib qualities: " << m_readPtr->getQualities());
        if ( m_readData.size() == 0 )
        {
            return;
        }

        int const startIndex = m_isForward ? 0 : m_readData.size() - 1;  // inclusive
        int const endIndex = m_isForward ? m_readData.size() - 1 : 0;    // inclusive
        int const direction = m_isForward ? 1 : -1;

        // find kickIndex — an index when error diploid kicks in first
        int kickIndex = startIndex;

        for ( /* empty */; /* empty */; kickIndex += direction )
        {
            if ( m_readData[kickIndex].errorPosterior() > 0.5 )
            {
                break;
            }
            if ( kickIndex == endIndex )  // no error state -- leave read alone
            {
                return;
            }
        }

        int readIndex =
            m_isForward ? m_readData[kickIndex].indexIntoRead() : m_readData[kickIndex].indexIntoRead() + kmerSize - 1;

        // ECHIDNA_LOG(SUPER_DEBUG, "Kick index: " << kickIndex);
        // ECHIDNA_LOG(SUPER_DEBUG, "Read index: " << readIndex);
        // ECHIDNA_LOG(SUPER_DEBUG, "Read data size: " << m_readData.size());
        do
        {
            // TODO(Edward, Hedvika) Is this safe from potential overflow?
            m_readPtr->qualities()[readIndex] = static_cast< char >( 2 );
            readIndex += direction;
        } while ( 0 <= readIndex && readIndex < m_readPtr->getLength() );

        // ECHIDNA_LOG(SUPER_DEBUG, "After recalib qualities: " << m_readPtr->getQualities());
    }
}
}
