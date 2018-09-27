// All content Copyright (C) 2018 Genomics plc
#include "mapping/hashMapper.hpp"
#include "utils/exceptions.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <string>
#include <map>
#include <numeric>

using namespace echidna::mapping;

int hashValueForSequence( const echidna::utils::BasePairSequence & sequence )
{
    const auto length = sequence.size();
    HashFunction hashFunction( sequence, length );
    return hashFunction.next( sequence.at( length - 1 ) );
}

BOOST_AUTO_TEST_CASE( test_hash_function_raises_for_kmer_of_length_1 )
{
    echidna::utils::BasePairSequence sequence( "A" );
    BOOST_CHECK_THROW( HashFunction( sequence, 1 ), echidna::utils::echidna_exception );
}

BOOST_AUTO_TEST_CASE( test_hash_function_raises_if_sequence_length_less_than_kmer_size )
{
    echidna::utils::BasePairSequence sequence( "AAAAAAAAAAA" );
    BOOST_CHECK_THROW( HashFunction( sequence, sequence.size() + 1 ), echidna::utils::echidna_exception );
}

BOOST_AUTO_TEST_CASE( test_hash_value_is_bijective_length_8 )
{
    const auto length = 8;
    const std::vector< char > letters = {'A', 'T', 'C', 'G'};
    const auto alphabetSize = letters.size();

    std::map< int, std::vector< std::string > > hashToKmers;

    const int numberOfKmers = int( pow( alphabetSize, length ) );
    for ( auto kmerIndex = 0; kmerIndex < numberOfKmers; ++kmerIndex )
    {
        std::string kmer( length, ' ' );
        for ( auto letterIndex = 0; letterIndex < length; ++letterIndex )
        {
            auto index = ( kmerIndex / ( int( pow( alphabetSize, letterIndex ) ) ) ) % alphabetSize;
            kmer[letterIndex] = letters[index];
        }

        const auto hashValueThisKmer = hashValueForSequence( kmer );

        for ( auto letter : letters )
        {
            const auto preAppendedKmer = std::string( 1, letter ) + kmer;
            HashFunction hashFunction( preAppendedKmer, length );
            hashFunction.next( preAppendedKmer[preAppendedKmer.size() - 2] );
            BOOST_REQUIRE_EQUAL( hashValueThisKmer, hashFunction.next( preAppendedKmer[preAppendedKmer.size() - 1] ) );
        }

        hashToKmers[hashValueThisKmer].push_back( kmer );
    }

    BOOST_CHECK_EQUAL( hashToKmers.size(), numberOfKmers );

    for ( const auto & pair : hashToKmers )
    {
        BOOST_REQUIRE_EQUAL( pair.second.size(), 1 );
    }
}

BOOST_AUTO_TEST_CASE( test_next_hash_value_repetitive_sequence )
{
    const auto kmerSize = 8;
    const echidna::utils::BasePairSequence sequence( "AAAAAAAA" );

    HashFunction hashFunction( sequence, kmerSize );

    auto firstHash = hashFunction.next( 'A' );
    BOOST_CHECK_EQUAL( firstHash, 0 );

    for ( int i = 0; i < kmerSize * 2; ++i )
    {
        auto nextHash = hashFunction.next( 'A' );
        BOOST_CHECK_EQUAL( firstHash, nextHash );
    }
}

BOOST_AUTO_TEST_CASE( test_maps_substring_back_to_correct_start_position_in_non_repetitive_sequence )
{
    const auto kmerSize = 8;
    const echidna::utils::BasePairSequence sequence( "ACGAAATTTAAAAACCCCCTCCGCCCAAAAACATGACCCAATTGG" );

    HashMapper hasher( sequence, kmerSize, 0 );

    for ( std::size_t start = 0; start < sequence.size() - kmerSize; ++start )
    {
        const echidna::utils::BasePairSequence read = sequence.substr( start, kmerSize );
        auto mappingPositions = hasher.mapSequence( read, boost::none );

        BOOST_REQUIRE_EQUAL( mappingPositions.size(), 1 );
        BOOST_CHECK_EQUAL( mappingPositions[0], start );
    }
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( test_fails_to_index_sequence_shorter_than_kmer_length )
{
    const auto kmerSize = 8;
    const echidna::utils::BasePairSequence sequence( kmerSize - 1, 'C' );

    BOOST_CHECK_THROW( HashMapper( sequence, kmerSize, 0 ), echidna::utils::echidna_exception );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( test_fails_to_map_read_shorter_than_kmer_length )
{
    const auto kmerSize = 8;
    const echidna::utils::BasePairSequence sequence( "ACGAAATTTAAAAACCCCCTCCGCCCAAAAACATGACCCAATTGG" );

    HashMapper hasher( sequence, kmerSize, 0 );

    const echidna::utils::BasePairSequence read = sequence.substr( 0, kmerSize - 1 );
    BOOST_CHECK_THROW( hasher.mapSequence( read, 0ul ), echidna::utils::echidna_exception );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( test_fails_to_index_sequence_longer_than_hash_table_size )
{
    const auto kmerSize = 8;
    const auto hashTableSize = std::pow( 4, kmerSize );
    const echidna::utils::BasePairSequence sequence( hashTableSize, 'A' );

    BOOST_CHECK_THROW( HashMapper( sequence, kmerSize, 0 ), echidna::utils::echidna_exception );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( test_fails_to_map_read_longer_than_hash_table_size )
{
    const auto kmerSize = 8;
    const auto hashTableSize = std::pow( 4, kmerSize );
    const echidna::utils::BasePairSequence sequence( 100, 'A' );

    HashMapper hasher( sequence, kmerSize, 0 );

    const echidna::utils::BasePairSequence read( hashTableSize + 1, 'A' );
    BOOST_CHECK_THROW( hasher.mapSequence( read, 0ul ), echidna::utils::echidna_exception );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( test_counts_kmer_matches_correctly_for_short_read_in_repetitive_sequence )
{
    const auto kmerSize = 8;
    const echidna::utils::BasePairSequence sequence( "AAAAAAAAAACCCCCCCCCCAAAAAAAAAA" );

    KmerMatches hasher( sequence, kmerSize, 0 );

    const echidna::utils::BasePairSequence read( kmerSize, 'A' );
    const auto kmerMatches = hasher.countKmerMatches( read );

    BOOST_CHECK_EQUAL( std::accumulate( kmerMatches.cbegin(), kmerMatches.cend(), 0 ), 6 );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( test_counts_kmer_matches_correctly_for_complete_read_of_repetitive_sequence )
{
    const auto kmerSize = 8;
    const echidna::utils::BasePairSequence sequence( "AAAAAAAAAA" );

    KmerMatches hasher( sequence, kmerSize, 0 );

    const echidna::utils::BasePairSequence read = sequence;
    auto counts = hasher.countKmerMatches( read );

    BOOST_CHECK_EQUAL( std::accumulate( counts.cbegin(), counts.cend(), 0 ), 3 );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( test_exact_matching_sequence_has_higher_matching_kmer_count_than_sequence_with_snp )
{
    const auto kmerSize = 8;
    const echidna::utils::BasePairSequence sequence( "AAAAAAAAAA" );
    const echidna::utils::BasePairSequence match( "AAAAAAAAAA" );
    const echidna::utils::BasePairSequence haplotypeWithSNP( "AAAATAAAAA" );

    KmerMatches hasher( sequence, kmerSize, 0 );

    const auto refMatches = hasher.countKmerMatches( match );
    const auto varMatches = hasher.countKmerMatches( haplotypeWithSNP );

    BOOST_CHECK_GT( *std::max_element( refMatches.cbegin(), refMatches.cend() ),
                    *std::max_element( varMatches.cbegin(), varMatches.cend() ) );

    BOOST_REQUIRE_EQUAL( refMatches.size(), 1 );
    BOOST_CHECK_EQUAL( refMatches[0], 3 );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(
    test_exact_matching_sequence_has_higher_matching_kmer_count_than_sequence_with_one_base_insertion )
{
    const auto kmerSize = 8;
    const echidna::utils::BasePairSequence sequence( "TTTTTTTTTAAAAAAAAAATTTTTTTTTTN" );
    const echidna::utils::BasePairSequence match( "TTTTTTTTTAAAAAAAAAATTTTTTTTTTN" );
    const echidna::utils::BasePairSequence haplotypeWithInsertion( "TTTTTTTTTAAAAAAAAAAATTTTTTTTTT" );

    KmerMatches hasher( sequence, kmerSize, 0 );

    const auto refMatches = hasher.countKmerMatches( match );
    const auto varMatches = hasher.countKmerMatches( haplotypeWithInsertion );

    BOOST_CHECK_GT( *std::max_element( refMatches.cbegin(), refMatches.cend() ),
                    *std::max_element( varMatches.cbegin(), varMatches.cend() ) );
}

BOOST_AUTO_TEST_CASE( test_should_only_report_one_match_if_read_and_sequence_the_same )
{
    const echidna::utils::BasePairSequence sequence( "AA" );
    const auto kmerSize = 2;
    KmerMatches matcher( sequence, kmerSize, 0 );

    const echidna::utils::BasePairSequence read( "AA" );
    const auto matches = matcher.countKmerMatches( read );

    const auto expectedNumMatches = sequence.size() - read.size() + 1;
    BOOST_REQUIRE_EQUAL( matches.size(), expectedNumMatches );
    BOOST_CHECK_EQUAL( matches[0], 1 );
}

BOOST_AUTO_TEST_CASE( test_should_report_no_matches_if_read_and_sequence_are_different )
{
    const echidna::utils::BasePairSequence sequence( "AA" );
    const auto kmerSize = 2;
    KmerMatches matcher( sequence, kmerSize, 0 );

    const echidna::utils::BasePairSequence read( "AT" );
    const auto matches = matcher.countKmerMatches( read );

    const auto expectedNumMatches = sequence.size() - read.size() + 1;
    BOOST_REQUIRE_EQUAL( matches.size(), expectedNumMatches );
    BOOST_CHECK_EQUAL( matches[0], 0 );
}

BOOST_AUTO_TEST_CASE( test_should_report_match_only_at_matching_index )
{
    const echidna::utils::BasePairSequence sequence( "AAT" );
    const auto kmerSize = 2;
    KmerMatches matcher( sequence, kmerSize, 0 );

    const echidna::utils::BasePairSequence read( "AT" );
    const auto matches = matcher.countKmerMatches( read );

    const auto expectedNumMatches = sequence.size() - read.size() + 1;
    BOOST_REQUIRE_EQUAL( matches.size(), expectedNumMatches );
    BOOST_CHECK_EQUAL( matches[0], 0 );
    BOOST_CHECK_EQUAL( matches[1], 1 );
}

BOOST_AUTO_TEST_CASE( test_should_report_number_of_matching_subkmers )
{
    const echidna::utils::BasePairSequence sequence( "AAAAT" );
    const auto kmerSize = 2;
    KmerMatches matcher( sequence, kmerSize, 0 );

    const echidna::utils::BasePairSequence read( "AAAA" );
    const auto matches = matcher.countKmerMatches( read );

    const auto expectedNumMatches = sequence.size() - read.size() + 1;
    BOOST_REQUIRE_EQUAL( matches.size(), expectedNumMatches );
    BOOST_CHECK_EQUAL( matches[0], 3 );
    BOOST_CHECK_EQUAL( matches[1], 2 );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( test_exact_matching_sequence_has_higher_matching_kmer_count_than_sequence_with_one_base_deletion )
{
    const auto kmerSize = 8;
    const echidna::utils::BasePairSequence sequence( "TTTTTTTTTAAAAAAAAAATTTTTTTTTT" );
    const echidna::utils::BasePairSequence match( "TTTTTTTTTAAAAAAAAAATTTTTTTTTT" );
    const echidna::utils::BasePairSequence haplotypeWithDeletion( "TTTTTTTTTAAAAAAAAATTTTTTTTTT" );

    KmerMatches hasher( sequence, kmerSize, 0 );

    const auto refMatches = hasher.countKmerMatches( match );
    const auto varMatches = hasher.countKmerMatches( haplotypeWithDeletion );

    BOOST_CHECK_GT( *std::max_element( refMatches.cbegin(), refMatches.cend() ),
                    *std::max_element( varMatches.cbegin(), varMatches.cend() ) );
}

//-------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( test_finds_multiple_mapping_positions_in_repetitive_sequence )
{
    const auto kmerSize = 8;
    //                                                0         1         2
    //                                                012345678901234567890123456789
    const echidna::utils::BasePairSequence sequence( "AAAAAAAAAACCCCCCCCCCAAAAAAAAAA" );

    HashMapper hasher( sequence, kmerSize, 0 );

    const echidna::utils::BasePairSequence read = sequence.substr( 0, kmerSize );
    auto mappingPositions = hasher.mapSequence( read, boost::none );

    BOOST_CHECK_EQUAL( mappingPositions.size(), 6 );
    BOOST_CHECK_EQUAL( mappingPositions[0], 22 );
    BOOST_CHECK_EQUAL( mappingPositions[1], 21 );
    BOOST_CHECK_EQUAL( mappingPositions[2], 20 );
    BOOST_CHECK_EQUAL( mappingPositions[3], 2 );
    BOOST_CHECK_EQUAL( mappingPositions[4], 1 );
    BOOST_CHECK_EQUAL( mappingPositions[5], 0 );
}

BOOST_AUTO_TEST_CASE( test_precondition_on_kmer_size )
{
    const echidna::utils::BasePairSequence sequence( "A" );
    const auto kmerSize = 1;
    BOOST_CHECK_THROW( HashMapper( sequence, kmerSize, 0 ), echidna::utils::echidna_exception );
}

BOOST_AUTO_TEST_CASE( test_sanity1 )
{
    const echidna::utils::BasePairSequence sequence( "AAA" );
    const auto kmerSize = 2;
    HashMapper hasher( sequence, kmerSize, 0 );

    const echidna::utils::BasePairSequence read( "AA" );
    auto mappingPositions = hasher.mapSequence( read, boost::none );

    BOOST_CHECK_EQUAL( mappingPositions.size(), 2 );
    BOOST_CHECK_EQUAL( mappingPositions[0], 1 );
    BOOST_CHECK_EQUAL( mappingPositions[1], 0 );
}

BOOST_AUTO_TEST_CASE( test_sanity2 )
{
    const echidna::utils::BasePairSequence sequence( "AAT" );
    const auto kmerSize = 2;
    HashMapper hasher( sequence, kmerSize, 0 );

    const echidna::utils::BasePairSequence read( "AA" );
    auto mappingPositions = hasher.mapSequence( read, boost::none );

    BOOST_CHECK_EQUAL( mappingPositions.size(), 1 );
    BOOST_CHECK_EQUAL( mappingPositions[0], 0 );
}

BOOST_AUTO_TEST_CASE( test_behaviour_when_no_decision_can_be_made )
{
    const echidna::utils::BasePairSequence sequence( "AAA" );
    const auto kmerSize = 2;
    HashMapper hasher( sequence, kmerSize, 0 );

    const echidna::utils::BasePairSequence read( "TT" );
    auto mappingPositions = hasher.mapSequence( read, boost::none );

    BOOST_REQUIRE_EQUAL( mappingPositions.size(), 0 );
}

BOOST_AUTO_TEST_CASE( test_sanity3 )
{
    const echidna::utils::BasePairSequence sequence( "ATAT" );
    const auto kmerSize = 2;
    HashMapper hasher( sequence, kmerSize, 0 );

    const echidna::utils::BasePairSequence read( "AT" );
    auto mappingPositions = hasher.mapSequence( read, boost::none );

    BOOST_CHECK_EQUAL( mappingPositions.size(), 2 );
    BOOST_CHECK_EQUAL( mappingPositions[0], 2 );
    BOOST_CHECK_EQUAL( mappingPositions[1], 0 );
}
