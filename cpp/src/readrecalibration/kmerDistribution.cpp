// All content Copyright (C) 2018 Genomics plc
#include "readrecalibration/kmerDistribution.hpp"
#include "alignment/cigarItems.hpp"
#include "caller/region.hpp"

namespace echidna
{
namespace corrector
{
    //-----------------------------------------------------------------------------------------

    KmerDistribution::KmerDistribution( const std::string & chromosomeLabel,
                                        const echidna::io::FastaFile & fa,
                                        const int readsStart,
                                        const int readsEnd )
        : m_firstReadStart( readsStart ), m_lastReadEnd( readsEnd )
    {
        const auto refSeqStart = m_firstReadStart - padding;
        const auto refSeqEnd = m_lastReadEnd + padding + kmerSize;

        const auto referenceSequenceString =
            fa.getSequence( caller::Region( chromosomeLabel, refSeqStart, refSeqEnd ) ).sequence();

        for ( auto pos = m_firstReadStart; pos < m_lastReadEnd; ++pos )
        {
            const auto paddedKmerStart = pos - refSeqStart - padding;
            const auto paddedKmerSize = kmerSize + 2 * padding;
            const auto extKmer = referenceSequenceString.substr( paddedKmerStart, paddedKmerSize );
            m_kmerDistribution.emplace_back( extKmer );
        }
    }

    //-----------------------------------------------------------------------------------------

    void KmerDistribution::resetErrorCountData( const double priorPerNucProbOfReadTurningIntoErrorState )
    {
        for ( auto & kmerDistrubtion : m_kmerDistribution )
        {
            kmerDistrubtion.resetErrorCountData( priorPerNucProbOfReadTurningIntoErrorState );
        }
    }

    //-----------------------------------------------------------------------------------------

    int KmerDistribution::posToIndex( const int pos ) const { return pos - m_firstReadStart; }

    //-----------------------------------------------------------------------------------------

    void KmerDistribution::updateErrorPosteriors()
    {
        for ( auto & kmerDistrubtion : m_kmerDistribution )
        {
            kmerDistrubtion.updateErrorProbabilities();
        }
    }

    //-----------------------------------------------------------------------------------------

    void KmerDistribution::accumulateErrorProbability( const int pos,
                                                       const double errorProbability,
                                                       const bool isForward )
    {
        const auto kmerIndex = this->posToIndex( pos );
        auto & kmer = m_kmerDistribution[kmerIndex];

        kmer.accumulateErrorProbability( errorProbability, isForward );
    }

    //-----------------------------------------------------------------------------------------

    const SiteKmerDistribution & KmerDistribution::getSiteKmerDistribution( const int pos ) const
    {
        const auto kmerIndex = this->posToIndex( pos );
        return m_kmerDistribution[kmerIndex];
    }

    //-----------------------------------------------------------------------------------------

    SiteKmerDistribution & KmerDistribution::getSiteKmerDistribution( const int pos )
    {
        const auto kmerIndex = this->posToIndex( pos );
        return m_kmerDistribution[kmerIndex];
    }

    //-----------------------------------------------------------------------------------------

    void KmerDistribution::updateKmerHistogram( const io::readPtr_t readPtr )
    {
        const auto referencePositions = readPtr->getReferencePositions();

        for ( std::size_t indexIntoRead = 0; indexIntoRead < referencePositions.size() - kmerSize + 1; ++indexIntoRead )
        {
            const int referencePosition = referencePositions[indexIntoRead];

            if ( referencePosition == alignment::emptyPos )
            {
                continue;
            }
            else
            {
                kmer_t< kmerSize > readKmer;
                {
                    const auto sequence = readPtr->sequence().substr( indexIntoRead, kmerSize );
                    std::copy( sequence.cbegin(), sequence.cend(), readKmer.data() );
                }

                SiteKmerDistribution & siteKmerDistribution = getSiteKmerDistribution( referencePosition );
                siteKmerDistribution.addKmer( readKmer );
            }
        }
    }

    //-----------------------------------------------------------------------------------------

    void KmerDistribution::finalise( const double probOfTrueKmerBeingRef )
    {
        for ( int pos = start(); pos < end(); ++pos )
        {
            std::unordered_map< kmer_t< kmerSize >, double, kmerhash_t< kmerSize > > kmerPriorMap;

            auto & siteKmerDistribution = m_kmerDistribution[posToIndex( pos )];
            const auto & refKmer = siteKmerDistribution.getReferenceKmer();
            double normalization = 0.0;

            for ( const auto & kmerCountPair : siteKmerDistribution.kmerCount() )
            {
                const auto & kmer = kmerCountPair.first;
                const auto & kmerCount = kmerCountPair.second;

                int priorCount = 0;

                if ( kmer == refKmer )
                {
                    priorCount = 1;
                }
                else
                {
                    for ( std::size_t i = 0; i < kmerSize; ++i )
                    {
                        priorCount -= ( kmer[i] == refKmer[i] ? 0 : 1 );  // subtract 1 for every mismatch
                    }
                    priorCount = std::max( -2, priorCount );  // Limit this number for sake of indels.
                }

                const double weight = std::max( 0, ( kmerCount + priorCount ) ) * ( kmerCount );

                if ( weight > 0.0 or kmer == refKmer )
                {
                    kmerPriorMap[kmer] = weight;
                    normalization += weight;
                }
            }

            const double refPrior = this->computeRefPrior( normalization, probOfTrueKmerBeingRef );

            kmerPriorMap[refKmer] += refPrior;
            normalization += refPrior;

            for ( auto & kmerPair : kmerPriorMap )
            {
                kmerPair.second /= normalization;

                if ( kmerPair.second > MINIMUM_KMER_PRIOR )
                {
                    siteKmerDistribution.addKmerPriorPair( kmerPair );
                }
            }
        }
    }

    //-----------------------------------------------------------------------------------------

    double KmerDistribution::computeRefPrior( const double normalisation, const double probOfTrueKmerBeingRef )
    {
        const double refPrior = normalisation / ( 1.0 - probOfTrueKmerBeingRef );

        if ( refPrior == 0.0 )
        {
            return 1.0;  // deal with the case of no (useful) data
        }
        else
        {
            return refPrior;
        }
    }

    //-----------------------------------------------------------------------------------------
}
}  // namespace echidna
