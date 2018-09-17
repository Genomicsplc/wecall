// All content Copyright (C) 2018 Genomics plc
#include "io/fastaFile.hpp"
#include "utils/exceptions.hpp"
#include "utils/referenceSequence.hpp"

#include <algorithm>
#include <stdexcept>
#include "caller/region.hpp"

namespace echidna
{
namespace io
{
    std::set< std::string > standardHumanChromosomes()
    {
        std::set< std::string > standardHumanChromosomes;
        for ( int chromIndex = 1; chromIndex <= 22; ++chromIndex )
        {
            standardHumanChromosomes.insert( std::to_string( chromIndex ) );
            standardHumanChromosomes.insert( "chr" + std::to_string( chromIndex ) );
        }
        std::set< std::string > nonAutosomes = {"X", "Y", "MT", "chrX", "chrY", "chrM"};
        standardHumanChromosomes.insert( nonAutosomes.begin(), nonAutosomes.end() );
        return standardHumanChromosomes;
    }

    std::string fastaIndexFileName( std::string fastaFileName ) { return fastaFileName + std::string( ".fai" ); }
    //-----------------------------------------------------------------------------------------
    // Implementation of FastaIndex functions.
    //-----------------------------------------------------------------------------------------

    //
    //               e.g., pass a flag --standardContigSet A, where A \el {human, tomato, …}
    //
    FastaIndex::FastaIndex( const std::string & fileName )
        : m_timer( std::make_shared< utils::Timer >( "IO", utils::fileMetaData( fileName ) ) )
    {
        // Constructor. Open index file, and parse contig information from
        // the file.
        const std::string indexFileName = std::string( fileName );
        ifstream theFile( indexFileName.c_str() );

        if ( not( theFile.is_open() ) )
        {
            throw utils::echidna_exception( "Could not open FASTA index file" );
        }

        parseContigInfoFromFile( theFile );
    }

    //-----------------------------------------------------------------------------------------

    void FastaIndex::parseContigInfoFromFile( ifstream & theFile )
    {
        utils::ScopedTimerTrigger scopedTimerTrigger( m_timer );
        std::string refName;
        int64_t seqLength = 0;
        int64_t start = 0;
        int64_t lineLength = 0;
        int64_t fullLineLength = 0;

        const auto standardHumanChroms = standardHumanChromosomes();

        while ( theFile >> refName >> seqLength >> start >> lineLength >> fullLineLength )
        {
            m_contigMap[refName] =
                std::make_shared< IndexTuple >( refName, seqLength, start, lineLength, fullLineLength );

            if ( standardHumanChroms.find( refName ) != standardHumanChroms.end() )
            {
                m_standardContigs.push_back( refName );
            }
        }
    }

    const itPtr_t FastaIndex::getIndexTuple( const std::string & refName ) const
    {
        // Retrieve the IndexTuple pointer corresponding to the specified
        // reference sequence contig.
        const auto it( m_contigMap.find( refName ) );
        ECHIDNA_ASSERT( it != m_contigMap.end(),
                        "Could not find refName: " + refName + " in FastaIndex::getIndexTuple function" );

        return it->second;
    }

    std::map< std::string, utils::Interval > FastaIndex::contigs() const
    {
        std::map< std::string, utils::Interval > contigs;
        for ( const auto & contigMapPair : m_contigMap )
        {
            const auto contigName = contigMapPair.first;
            contigs.insert( std::make_pair( contigName, utils::Interval( 0L, this->getContigLength( contigName ) ) ) );
        }
        return contigs;
    }

    //-----------------------------------------------------------------------------------------
    // Implementation of FastaFile functions.
    //-----------------------------------------------------------------------------------------

    FastaFile::FastaFile( const std::string & fileName )
        : m_cacheSequence( nullptr ),
          m_file( fileName.c_str() ),
          m_indexFile( fastaIndexFileName( fileName ) ),
          m_timer( std::make_shared< utils::Timer >( "IO", utils::fileMetaData( fileName ) ) )
    {
        if ( not( m_file.is_open() ) )
        {
            throw utils::echidna_exception( "Could not open FASTA file" );
        }
    }

    //-----------------------------------------------------------------------------------------

    std::pair< int64_t, int64_t > FastaFile::computeStartAndEndPosOfSequenceInFile(
        const caller::Region & region ) const
    {

        const auto theTuple = m_indexFile.getIndexTuple( region.contig() );

        const auto contigInterval = utils::Interval( 0L, this->m_indexFile.getContigLength( region.contig() ) );

        ECHIDNA_ASSERT( contigInterval.contains( region.interval() ),
                        "Region " + region.toString() + " is not valid as not contained contig" );

        const auto nEndOfLineBytesBeforeStart =
            ( ( theTuple->m_fullLineLength - theTuple->m_lineLength ) * region.start() ) / theTuple->m_lineLength;
        const auto nEndOfLineBytesBeforeEnd =
            ( ( theTuple->m_fullLineLength - theTuple->m_lineLength ) * region.end() ) / theTuple->m_lineLength;

        const auto seqStartPosInFile =
            theTuple->m_startPos + std::max( int64_t( 0 ), region.start() ) + nEndOfLineBytesBeforeStart;
        const auto seqEndPosInFile =
            theTuple->m_startPos + std::min( theTuple->m_seqLength, region.end() ) + nEndOfLineBytesBeforeEnd;

        return std::make_pair( seqStartPosInFile, seqEndPosInFile );
    }

    //-----------------------------------------------------------------------------------------

    void FastaFile::getPaddedSequenceFromFile( std::string * seq, const caller::Region & region ) const
    {
        auto const contigLength = m_indexFile.getContigLength( region.contig() );
        const auto numPaddedLeft = std::max( -region.start(), 0L );
        const auto numPaddedRight = std::max( region.end() - contigLength, 0L );
        const auto unpaddedStartPos = region.start() + numPaddedLeft;
        const auto unpaddedEndPos = region.end() - numPaddedRight;

        seq->append( int64_to_sizet( numPaddedLeft ), constants::gapChar );
        const caller::Region unpaddedRegion( region.contig(), unpaddedStartPos, unpaddedEndPos );

        {
            utils::ScopedTimerTrigger scopedTimerTrigger( m_timer );
            // Retrieve the specified sequence from the FASTA file, and cache it to a
            // string buffer.
            const auto fileStartEnd = computeStartAndEndPosOfSequenceInFile( unpaddedRegion );

            const int64_t bufferSize = fileStartEnd.second - fileStartEnd.first;
            auto tempBuffer = std::unique_ptr< char[] >( new char[bufferSize] );

            m_file.seekg( fileStartEnd.first );
            m_file.read( tempBuffer.get(), bufferSize );

            for ( int64_t i = 0; i < bufferSize; ++i )
            {
                const char theChar = tempBuffer[i];

                if ( theChar != '\n' )
                {
                    seq->push_back( std::toupper( theChar ) );
                }
            }
        }
        seq->append( int64_to_sizet( numPaddedRight ), constants::gapChar );
    }

    //-----------------------------------------------------------------------------------------

    void FastaFile::cacheSequence( const caller::Region & region )
    {
        std::string cache;
        this->getPaddedSequenceFromFile( &cache, region );

        m_cacheSequence = std::make_shared< utils::ReferenceSequence >( region, cache );
    }

    utils::ReferenceSequence FastaFile::getSequence( const caller::Region & region ) const
    {
        if ( m_cacheSequence != nullptr and m_cacheSequence->region().contains( region ) )
        {
            return m_cacheSequence->subseq( region );
        }
        else
        {
            std::string seq;
            this->getPaddedSequenceFromFile( &seq, region );
            return utils::ReferenceSequence( region, seq );
        }
    }

    //-----------------------------------------------------------------------------------------
}
}
