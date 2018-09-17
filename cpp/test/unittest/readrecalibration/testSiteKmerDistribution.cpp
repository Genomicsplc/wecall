#include "io/fastaFile.hpp"
#include "readrecalibration/siteKmerDistribution.hpp"
#include "readrecalibration/commonTypes.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( shouldConstructReferenceStringFromPaddedReference )
{
    // Given
    //                             123456789
    std::string paddedReference = "ATCGCTCTG";
    BOOST_CHECK_EQUAL( paddedReference.size(), 2 * echidna::corrector::padding + echidna::corrector::kmerSize );

    // When
    echidna::corrector::SiteKmerDistribution siteKmerDistribution( paddedReference );

    // Then
    BOOST_CHECK_EQUAL( echidna::corrector::show_string( siteKmerDistribution.getReferenceKmer() ),
                       paddedReference.substr( 1, paddedReference.size() - 2 ) );
    BOOST_CHECK_EQUAL( echidna::corrector::show_string( siteKmerDistribution.getExtReferenceKmer() ), paddedReference );
}

BOOST_AUTO_TEST_CASE( shouldComputePErrorAsPriorOnBeingReset )
{
    std::string paddedReference = "ATCGCTCTG";
    echidna::corrector::SiteKmerDistribution siteKmerDistribution( paddedReference );
    auto priorPerNucProbOfReadTurningIntoErrorState = 0.661238;

    siteKmerDistribution.resetErrorCountData( priorPerNucProbOfReadTurningIntoErrorState );
    siteKmerDistribution.updateErrorProbabilities();

    BOOST_CHECK_CLOSE( siteKmerDistribution.pError( true ), priorPerNucProbOfReadTurningIntoErrorState, 1e-8 );
    BOOST_CHECK_CLOSE( siteKmerDistribution.pError( false ), priorPerNucProbOfReadTurningIntoErrorState, 1e-8 );
}
