// All content Copyright (C) 2018 Genomics plc
#include "readrecalibration/siteKmerDistribution.hpp"
#include <cassert>

namespace wecall
{
namespace corrector
{
    SiteKmerDistribution::SiteKmerDistribution( wecall::utils::BasePairSequence paddedKmer )
    {
        assert( paddedKmer.size() == kmerSize + 2 * padding );

        std::copy( paddedKmer.cbegin() + padding, paddedKmer.cend() - padding, m_referenceKmer.data() );
        std::copy( paddedKmer.cbegin(), paddedKmer.cend(), m_extReferenceKmer.data() );
    }

    void SiteKmerDistribution::resetErrorCountData( double priorPerNucProbOfReadTurningIntoErrorState )
    {
        double const pseudoCount = 1;  // modestly conservative setPrior: pretend to see one well-behaving read
        m_errorCountData[0] = {pseudoCount, pseudoCount * priorPerNucProbOfReadTurningIntoErrorState};
        m_errorCountData[1] = {pseudoCount, pseudoCount * priorPerNucProbOfReadTurningIntoErrorState};
    }

    void SiteKmerDistribution::accumulateErrorProbability( double errorProbability, bool isForward )
    {
        int fwbwIndex = isForward ? 0 : 1;

        assert( errorProbability >= -1.0e-10 );
        assert( m_errorCountData[fwbwIndex].errorOpportunity > 0.0 );
        assert( m_errorCountData[fwbwIndex].errorCount > 0.0 );

        m_errorCountData[fwbwIndex].errorOpportunity += 1;
        m_errorCountData[fwbwIndex].errorCount += errorProbability;
    }

    void SiteKmerDistribution::updateErrorProbabilities()
    {
        m_pErrorForward = m_errorCountData[0].errorCount / m_errorCountData[0].errorOpportunity;
        m_pErrorBackward = m_errorCountData[1].errorCount / m_errorCountData[1].errorOpportunity;
    }
}
}
