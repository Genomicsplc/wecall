// All content Copyright (C) 2018 Genomics plc
#ifndef SITE_KMER_DISTRIBUTION_HPP
#define SITE_KMER_DISTRIBUTION_HPP

#include <string>
#include <map>
#include <vector>
#include <unordered_map>
#include <utility>
#include <utils/sequence.hpp>

#include "readrecalibration/errorCorrectionParameters.hpp"
#include "readrecalibration/commonTypes.hpp"

namespace echidna
{
namespace corrector
{
    class SiteKmerDistribution
    {

    public:
        SiteKmerDistribution( echidna::utils::BasePairSequence paddedKmer );

        std::size_t size() const { return m_kmerCount.size(); }

        double pError( bool isForward ) const { return isForward ? m_pErrorForward : m_pErrorBackward; }

        const kmer_t< kmerSize > & getReferenceKmer() const { return m_referenceKmer; }

        const extKmer_t< kmerSize, padding > & getExtReferenceKmer() const { return m_extReferenceKmer; }

        const std::unordered_map< kmer_t< kmerSize >, int, kmerhash_t< kmerSize > > & kmerCount() const
        {
            return m_kmerCount;
        }

        const std::vector< std::pair< kmer_t< kmerSize >, double > > & kmerPrior() const { return m_kmerPrior; }

        void addKmerPriorPair( std::pair< kmer_t< kmerSize >, double > kmerPriorPair )
        {
            m_kmerPrior.push_back( kmerPriorPair );
        }

    public:
        void addKmer( const kmer_t< kmerSize > & kmer )
        {
            auto & countThisKmer = m_kmerCount[kmer];
            ++countThisKmer;
        }

        void resetErrorCountData( double priorPerNucProbOfReadTurningIntoErrorState );

        void accumulateErrorProbability( double errorProbability, bool isForward );

        void updateErrorProbabilities();

    private:
        // TODO(semen): consider collapsing m_referenceKmer and m_extReferenceKmer
        kmer_t< kmerSize > m_referenceKmer;                 ///< kmer in the genome reference at this location
        extKmer_t< kmerSize, padding > m_extReferenceKmer;  ///< padded kmer in genome reference at this location

        std::unordered_map< kmer_t< kmerSize >, int, kmerhash_t< kmerSize > >
            m_kmerCount;  ///< histogram of kmers observed at this location
        std::vector< std::pair< kmer_t< kmerSize >, double > > m_kmerPrior;

        double m_pErrorForward = 0;   ///< probability of moving into the error state on forward reads
        double m_pErrorBackward = 0;  ///< probability of moving into the error state on backward reads
        std::vector< ErrorCountData > m_errorCountData = {{0, 0}, {0, 0}};
    };
}  // namespace corrector
}  // namespace echidna

#endif
