#include "stats/functions.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "stats/models.hpp"
#include "common.hpp"

BOOST_AUTO_TEST_CASE( testGypergeometricFunctionThreeFTwo )
{
    const int alpha = 20;
    const int beta = 20;

    int k = 0;
    int n = 30;
    BOOST_CHECK_CLOSE( echidna::stats::threeFtwo( alpha + k + 1, -n + k + 1, k + 2, -beta - n + k + 2 ),
                       1013883283.0 / 7400.0, 1e-3 );

    n = 100;
    BOOST_CHECK_CLOSE( echidna::stats::threeFtwo( alpha + k + 1, -n + k + 1, k + 2, -beta - n + k + 2 ),
                       17075323167913883149.0 / 27898000.0, 1e-3 );

    n = 150;
    BOOST_CHECK_CLOSE( echidna::stats::threeFtwo( alpha + k + 1, -n + k + 1, k + 2, -beta - n + k + 2 ),
                       2405633174463118093.0 / 6000.0, 1e-3 );

    n = 200;
    BOOST_CHECK_CLOSE( echidna::stats::threeFtwo( alpha + k + 1, -n + k + 1, k + 2, -beta - n + k + 2 ),
                       9915559267243212161923.0 / 186000.0, 1e-3 );

    k = 1;
    n = 30;
    BOOST_CHECK_CLOSE( echidna::stats::threeFtwo( alpha + k + 1, -n + k + 1, k + 2, -beta - n + k + 2 ),
                       4055503532.0 / 187775.0, 1e-3 );
}

BOOST_AUTO_TEST_CASE( test_beta_binomial_example )
{
    const int alpha = 180;
    const int beta = 20;
    const int totalCoverage = 3;

    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDF( 0, totalCoverage, alpha, beta ), 0.00113787, 1e-3 );
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDF( 1, totalCoverage, alpha, beta ), 0.0290675, 1e-3 );
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDF( 2, totalCoverage, alpha, beta ), 0.269795, 1e-3 );
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDF( 3, totalCoverage, alpha, beta ), 1.0, 1e-3 );
}

BOOST_AUTO_TEST_CASE( testBetaBinomialForExtremeCasesInReferenceCalling )
{
    // alpha and beta as used for reference calling
    const int alpha = 20;
    const int beta = 20;

    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDF( 0, 30, alpha, beta ), 259.0 / 434521666.0, 1e-3 );
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDF( 0, 80, alpha, beta ), 2.8072527061377133e-12, 7 );

    // machine precision lower bound is used
    // in reality this should be 9.72123×10^-14
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDF( 0, 100, alpha, beta ), 6.5114324144465035e-14, 1e-3 );
    // in reality this should be 5.47187×10^-15
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDF( 0, 120, alpha, beta ), 6.5478967258053169e-14, 1e-3 );
    // in reality this should be 2.83125×10^-15
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDF( 0, 125, alpha, beta ), 6.5560611247093511e-14, 1e-3 );
}

BOOST_AUTO_TEST_CASE( testBetaBinomialStableForExtremeCasesInReferenceCalling )
{
    // alpha and beta as used for reference calling
    const int alpha = 20;
    const int beta = 20;

    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDFStable( 0, 30, alpha, beta ), 5.9605779867411002e-07, 1e-3 );
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDFStable( 0, 80, alpha, beta ), 2.8072527061377133e-12, 0.1 );
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDFStable( 0, 100, alpha, beta ), 9.72123e-14, 1 );
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDFStable( 0, 125, alpha, beta ), 2.83125e-15, 17 );

    // approximation becomes very bad so return a fixed value
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDFStable( 0, 150, alpha, beta ), 1e-15, 1e-3 );
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDFStable( 0, 175, alpha, beta ), 1e-15, 1e-3 );
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDFStable( 0, 200, alpha, beta ), 1e-15, 1e-3 );
}

BOOST_AUTO_TEST_CASE( testBetaBinomialForReferenceCalls )
{
    // alpha and beta as used for reference calling
    const int alpha = 20;

    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDFForReferenceCalls( 20, alpha ), 2.4663342453545845e-05, 1e-3 );
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDFForReferenceCalls( 40, alpha ), 2.5994125629758003e-08, 1e-3 );
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDFForReferenceCalls( 60, alpha ), 1.6074021186066192e-10, 1e-3 );
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDFForReferenceCalls( 80, alpha ), 2.8072527061377133e-12, 0.1 );

    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDFForReferenceCalls( 100, alpha ), 9.72123e-14, 1.1 );
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDFForReferenceCalls( 125, alpha ), 2.83125e-15, 2.1 );

    // approximation becomes very bad so return a fixed value
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDFForReferenceCalls( 150, alpha ), 1e-15, 1e-3 );
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDFForReferenceCalls( 175, alpha ), 1e-15, 1e-3 );
    BOOST_CHECK_CLOSE( echidna::stats::betaBinomialCDFForReferenceCalls( 200, alpha ), 1e-15, 1e-3 );
}

BOOST_AUTO_TEST_CASE( shouldComputeRootMeanSquare )
{
    std::vector< double > input = {1.0, 2.0, 3.0};
    const auto expectedResult = 2.160246899469287;
    BOOST_CHECK_CLOSE( expectedResult, echidna::stats::rootMeanSquare( input ), 1e-5 );
}

BOOST_AUTO_TEST_CASE( shouldReturnZeroIfVariantSupportIsZero )
{
    BOOST_CHECK( std::isnan( echidna::stats::variantSupportPerRead( 1.0, 1.0, 0 ) ) );
}

BOOST_AUTO_TEST_CASE( shouldReturnCorrectRatioIfVariantSupportNonZero )
{
    const auto expectedResult = 1.0;
    BOOST_CHECK_CLOSE( expectedResult, echidna::stats::variantSupportPerRead( 1.0, 1.0, 2 ), 1e-5 );
}

BOOST_AUTO_TEST_CASE( shouldReturnNaNIfPosteriorAndPriorTooLarge )
{
    const auto result = echidna::stats::variantSupportPerRead( 1.0, constants::maxPhredScore - 1.0, 2 );
    BOOST_CHECK( std::isnan( result ) );
}

BOOST_AUTO_TEST_CASE( shouldReturnNaNIfPosteriorTooLarge )
{
    const auto result = echidna::stats::variantSupportPerRead( 0.0, constants::maxPhredScore, 2 );
    BOOST_CHECK( std::isnan( result ) );
}

BOOST_AUTO_TEST_CASE( shouldReturnNaNIfPriorTooLarge )
{
    const auto result = echidna::stats::variantSupportPerRead( constants::maxPhredScore, 0.0, 2 );
    BOOST_CHECK( std::isnan( result ) );
}

BOOST_AUTO_TEST_CASE( shouldReturnNaNIfCombinedSumOfPriorAndPosteriorTooLarge )
{
    const auto result = echidna::stats::variantSupportPerRead( constants::maxPhredScore, 0.0, 2 );
    BOOST_CHECK( std::isnan( result ) );
}

BOOST_AUTO_TEST_CASE( shouldClipToMaxPhredScoreIfProbErrIsTooLow )
{
    const auto lowProbErr = 1e-301;
    BOOST_CHECK_EQUAL( constants::maxPhredScore, echidna::stats::toPhredQ( lowProbErr ) );
}

BOOST_AUTO_TEST_CASE( shouldReturnComputedPhredScoreForLowProbability )
{
    const auto lowProbErr = 1e-199;
    BOOST_CHECK_EQUAL( 1990, echidna::stats::toPhredQ( lowProbErr ) );
}

BOOST_AUTO_TEST_CASE( shouldReturnMaxPhredScoreIfProbErrIsZero )
{
    const auto lowProbErr = 0.0;
    BOOST_CHECK_EQUAL( constants::maxPhredScore, echidna::stats::toPhredQ( lowProbErr ) );
}

BOOST_AUTO_TEST_CASE( shouldReturnMaxPhredScoreIfProbErrIsOne )
{
    const auto lowProbErr = 1.0;
    BOOST_CHECK_EQUAL( 0.0, echidna::stats::toPhredQ( lowProbErr ) );
}

BOOST_AUTO_TEST_CASE( shouldReturnMaxPhredScoreIfProbErrIsOneTenth )
{
    const auto lowProbErr = 0.1;
    BOOST_CHECK_EQUAL( 10.0, echidna::stats::toPhredQ( lowProbErr ) );
}

BOOST_AUTO_TEST_CASE( shouldReturnProbOneTenthFromPhredScoreTen )
{
    const auto phredScore = 10.0;
    BOOST_CHECK_EQUAL( 0.1, echidna::stats::fromPhredQ( phredScore ) );
}

BOOST_AUTO_TEST_CASE( shouldReturnProbOneFromPhredScoreOfZero )
{
    const auto phredScore = 0.0;
    BOOST_CHECK_EQUAL( 1.0, echidna::stats::fromPhredQ( phredScore ) );
}

BOOST_AUTO_TEST_CASE( shouldReturnLowProbFromHighPhredScore )
{
    const auto phredScore = 1000.0;
    BOOST_CHECK_EQUAL( 1e-100, echidna::stats::fromPhredQ( phredScore ) );
}
