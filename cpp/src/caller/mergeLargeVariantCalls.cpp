#include "mergeLargeVariantCalls.hpp"

namespace echidna
{
namespace caller
{
    callVector_t lvcMerger::removeReferenceCalls( const callVector_t & calls,
                                                  const std::vector< std::size_t > & thisAreaPloidies,
                                                  const std::vector< std::size_t > & defaultPloidies ) const
    {
        callVector_t filteredRefCalls;
        for ( auto & call : calls )
        {
            if ( not call.isRefCall() )
            {
                filteredRefCalls.push_back( call );
            }
            else
            {
                // keep reference calls if a sample does not have the large variant
                bool keepRefCall = true;
                for ( size_t sampleIndex = 0u; sampleIndex < defaultPloidies.size(); sampleIndex++ )
                    if ( thisAreaPloidies[sampleIndex] < defaultPloidies[sampleIndex] )
                    {
                        keepRefCall = false;
                        continue;
                    }
                if ( keepRefCall )
                {
                    filteredRefCalls.push_back( call );
                }
            }
        }
        return filteredRefCalls;
    }

    std::vector< size_t > lvcMerger::getOverlappingLvcCallIndices( const Call & call,
                                                                   const callVector_t & lvcCalls ) const
    {
        std::vector< size_t > overlappingCallIndices = {};

        for ( size_t callIndex = 0; callIndex < lvcCalls.size(); callIndex++ )
        {
            if ( not lvcCalls[callIndex].isRefCall() and
                 lvcCalls[callIndex].var->region().overlaps( call.var->region() ) )
            {
                overlappingCallIndices.push_back( callIndex );
            }
        }
        return overlappingCallIndices;
    }

    // As the large variant cluster and the other variant cluster where called seperately, the genotype is out of sync.
    // This function will adapt the genotype of the calls in the other variant cluster if they overlap the large
    // variant.
    void lvcMerger::mergeAndCorrectGenotypes( callVector_t & otherCalls,
                                              callVector_t & lvcCalls,
                                              const std::vector< std::size_t > & thisAreaPloidies,
                                              const std::vector< std::size_t > & defaultPloidies ) const
    {
        for ( auto & call : otherCalls )
        {
            for ( std::size_t sampleIndex = 0; sampleIndex < defaultPloidies.size(); ++sampleIndex )
            {
                if ( call.samples[sampleIndex].genotypeCalls.size() == 2 )
                {
                    assert( thisAreaPloidies[sampleIndex] == 2 );
                    continue;
                }

                if ( call.samples[sampleIndex].genotypeCalls.size() == 0 )
                {
                    auto overlappingLvcCallIndices = this->getOverlappingLvcCallIndices( call, lvcCalls );
                    ECHIDNA_ASSERT( overlappingLvcCallIndices.size() >= 1,
                                    "No overlapping call found for " + call.var->toString() );

                    ECHIDNA_ASSERT( thisAreaPloidies[sampleIndex] == 0, "wrong ploidy" );
                    call.samples[sampleIndex].genotypeCalls = genoCall_t( defaultPloidies[sampleIndex], Call::UNKNOWN );
                    continue;
                }

                if ( call.samples[sampleIndex].genotypeCalls.size() == 1 )
                {
                    ECHIDNA_ASSERT( thisAreaPloidies[sampleIndex] == 1, "wrong ploidy" );

                    auto overlappingLvcCallIndices = this->getOverlappingLvcCallIndices( call, lvcCalls );
                    ECHIDNA_ASSERT( overlappingLvcCallIndices.size() >= 1,
                                    "No overlapping call found for " + call.var->toString() );

                    for ( size_t & lvIndex : overlappingLvcCallIndices )
                    {

                        if ( lvcCalls[lvIndex].samples[sampleIndex].genotypeCalls[0] == Call::VAR )
                        {
                            ECHIDNA_ASSERT(
                                lvcCalls[lvIndex].samples[sampleIndex].genotypeCalls[1] == Call::REF ||
                                    lvcCalls[lvIndex].samples[sampleIndex].genotypeCalls[1] == Call::UNKNOWN,
                                "Unexpected genotype found!" +
                                    std::to_string( lvcCalls[lvIndex].samples[sampleIndex].genotypeCalls[1] ) );
                            call.samples[sampleIndex].genotypeCalls.insert(
                                call.samples[sampleIndex].genotypeCalls.begin(), Call::UNKNOWN );
                            lvcCalls[lvIndex].samples[sampleIndex].genotypeCalls[1] = Call::UNKNOWN;
                        }
                        else if ( lvcCalls[lvIndex].samples[sampleIndex].genotypeCalls[1] == Call::VAR )
                        {
                            ECHIDNA_ASSERT(
                                lvcCalls[lvIndex].samples[sampleIndex].genotypeCalls[0] == Call::REF ||
                                    lvcCalls[lvIndex].samples[sampleIndex].genotypeCalls[0] == Call::UNKNOWN,
                                "Unexpected genotype found!" +
                                    std::to_string( lvcCalls[lvIndex].samples[sampleIndex].genotypeCalls[0] ) );
                            call.samples[sampleIndex].genotypeCalls.push_back( Call::UNKNOWN );
                            lvcCalls[lvIndex].samples[sampleIndex].genotypeCalls[0] = caller::Call::UNKNOWN;
                        }
                        else
                        {
                            // it can happen that the large variant is not called for a sample when there is a second
                            // large variant - do nothing in this case
                        }
                    }
                }
                else
                {
                    ECHIDNA_ASSERT( false, "Unexpected structure of large variant call found." );
                }
            }
        }
    }
}
}
