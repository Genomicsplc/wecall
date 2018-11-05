// All content Copyright (C) 2018 Genomics plc
#ifndef KMER_DISTRIBUTION_HPP
#define KMER_DISTRIBUTION_HPP

#include <string>
#include <map>
#include <utility>

#include "io/read.hpp"
#include "io/fastaFile.hpp"

#include "readrecalibration/siteKmerDistribution.hpp"

namespace wecall
{
namespace corrector
{
    class KmerDistribution
    {

    public:
        KmerDistribution( const std::string & chromosomeLabel,
                          const io::FastaFile & fa,
                          const int readsStart,
                          const int readsEnd );

        int start() const { return m_firstReadStart; }
        int end() const { return m_lastReadEnd; }

        void updateKmerHistogram( const io::readPtr_t readPtr );
        void finalise( const double probOfTrueKmerBeingRef );
        void resetErrorCountData( const double priorPerNucProbOfReadTurningIntoErrorState );
        void updateErrorPosteriors();

        void accumulateErrorProbability( const int pos, const double errorProbability, const bool isForward );

        const SiteKmerDistribution & getSiteKmerDistribution( const int pos ) const;
        SiteKmerDistribution & getSiteKmerDistribution( const int pos );

    private:
        double computeRefPrior( const double normalisation, const double probOfTrueKmerBeingRef );

        int posToIndex( const int pos ) const;

    private:
        const int m_firstReadStart;
        const int m_lastReadEnd;

        std::vector< SiteKmerDistribution > m_kmerDistribution;
    };
}  // namespace corrector
}  // namespace wecall

#endif
