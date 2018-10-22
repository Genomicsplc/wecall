// All content Copyright (C) 2018 Genomics plc
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "caller/haplotypeLikelihoods.hpp"
#include "alignment/aligner.hpp"
#include "utils/matrix.hpp"

namespace wecall
{
namespace caller
{
    std::vector< double > computeHaplotypeFrequencies( const utils::matrix_t & haplotypeLikelihoods )
    {
        std::vector< double > haplotypeFrequencies( haplotypeLikelihoods.size2(), 0.0 );
        for ( std::size_t haplotypeIndex = 0; haplotypeIndex != haplotypeLikelihoods.size2(); ++haplotypeIndex )
        {
            const auto column = utils::matrixColumn_t( haplotypeLikelihoods, haplotypeIndex );
            const auto hapFreq = std::accumulate( column.begin(), column.end(), 0.0 );
            haplotypeFrequencies[haplotypeIndex] = hapFreq;
        }
        return haplotypeFrequencies;
    }

    utils::matrix_t computeHaplotypeLikelihoods( const variant::HaplotypeVector & haplotypes,
                                                 const io::RegionsReads & readRange )
    {
        const auto nReads = static_cast< std::size_t >( std::distance( readRange.begin(), readRange.end() ) );
        const auto haplotypeSeqStart = haplotypes.paddedReferenceSequence()->start();

        utils::matrix_t readScores( nReads, haplotypes.size() );

        for ( std::size_t haplotypeIndex = 0; haplotypeIndex < haplotypes.size(); ++haplotypeIndex )
        {
            const auto & paddedHaplotypeSequences = haplotypes[haplotypeIndex].paddedSequences();

            std::vector< mapping::HashMapper > hashMappers;
            std::vector< alignment::GAlign > aligners;
            for ( const auto & paddedHaplotypeSequence : paddedHaplotypeSequences )
            {
                hashMappers.emplace_back( paddedHaplotypeSequence, constants::needlemanWunschPadding,
                                          constants::needlemanWunschPadding );
                aligners.emplace_back(
                    paddedHaplotypeSequence, constants::gapExtendPenalty, constants::nucleotidePrior,
                    alignment::computeGapOpen( paddedHaplotypeSequence, errorModels::illuminaErrorModel ) );
            }

            std::size_t readIndex = 0;
            for ( const auto & read : readRange )
            {
                const auto hintPosition = read.getStartPos() - haplotypeSeqStart;
                double maxReadScore = 0.0;
                for ( std::size_t sequenceIndex = 0; sequenceIndex < paddedHaplotypeSequences.size(); ++sequenceIndex )
                {
                    maxReadScore = std::max(
                        maxReadScore, computeLikelihoodForReadAndHaplotype(
                                          read, hintPosition, hashMappers[sequenceIndex], aligners[sequenceIndex] ) );
                }

                readScores( readIndex, haplotypeIndex ) = maxReadScore;
                ++readIndex;
            }
        }

        int64_t maxMappingQual = 0L;
        for ( const auto & read : readRange )
        {
            maxMappingQual = std::max( maxMappingQual, read.getMappingQuality() );
        }

        const double maxDifference = stats::fromPhredQ( maxMappingQual );

        utils::smoothLowOutliers( readScores, maxDifference );

        if ( false )
        {
            const auto debugFunc = [&haplotypes, &readScores, &readRange]() -> std::string
            {
                std::stringstream debugMessage;
                debugMessage << "Haplotypes order:- ";
                for ( const auto & haplotype : haplotypes )
                {
                    debugMessage << std::endl
                                 << haplotype.toString();
                }

                debugMessage << std::endl;
                for ( std::size_t hap_index = 0; hap_index < readScores.size2(); ++hap_index )
                {
                    debugMessage << "\tHaplotype scores for : " << haplotypes[hap_index].toString() << std::endl;
                    for ( std::size_t read_index = 0; read_index < readScores.size1(); ++read_index )
                    {
                        debugMessage << readScores( read_index, hap_index ) << ", ";
                    }
                    debugMessage << std::endl;
                }
                return debugMessage.str();
            };
            ECHIDNA_LOG( SUPER_DEBUG, debugFunc() );
        }

        return readScores;
    }
}
}
