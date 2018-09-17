// All content Copyright (C) 2018 Genomics plc
#include <cassert>
#include "io/bamFileIterator.hpp"

namespace echidna
{
namespace io
{
    AbstractBamFileIterator::AbstractBamFileIterator( bam_fetch_iterator_t * bamIterator,
                                                      utils::referenceSequencePtr_t refSequence,
                                                      utils::timerPtr_t timer )
        : m_refSequence( refSequence ), m_bamIterator( bamIterator ), m_timer( timer )
    {
        m_bamRecordPtr = bam_fetch_iterate( m_bamIterator );
    }

    AbstractBamFileIterator::~AbstractBamFileIterator()
    {
        bam_cleanup_fetch_iterator( m_bamIterator );
        free( m_bamIterator );
    }

    void AbstractBamFileIterator::next() { m_bamRecordPtr = bam_fetch_iterate( m_bamIterator ); }

    bool AbstractBamFileIterator::hasReadData() { return ( m_bamRecordPtr != nullptr ); }

    //-----------------------------------------------------------------------------------------

    BamFileIterator::BamFileIterator( bam_fetch_iterator_t * bamIterator,
                                      std::map< std::string, std::string > samplesByID,
                                      utils::referenceSequencePtr_t refSequence,
                                      utils::timerPtr_t timer )
        : AbstractBamFileIterator( bamIterator, refSequence, timer ), m_samplesByID( samplesByID )
    {
    }

    BamFileIterator::~BamFileIterator() {}

    std::pair< std::string, readPtr_t > BamFileIterator::getReadData()
    {
        if ( not this->hasReadData() )
        {
            throw utils::echidna_exception( "Tried to get read data without checking if there were more to get" );
        }

        const uint8_t * rgID = bam_aux_get( m_bamRecordPtr, "RG" );
        assert( rgID );

        const char * rgID_char_ptr = bam_aux2Z( rgID );
        assert( rgID_char_ptr );
        const std::string rgIDStr = rgID_char_ptr;

        auto it = m_samplesByID.find( rgIDStr );

        if ( it == m_samplesByID.end() )
        {
            throw utils::echidna_exception( "Found read without matching sample name in BAM header!" );
        }

        std::string sampleName = it->second;

        readPtr_t read = std::make_shared< Read >( m_bamRecordPtr, m_refSequence );

        return std::make_pair( sampleName, read );
    }

    //-----------------------------------------------------------------------------------------

    BamFileWithoutReadGroupIterator::BamFileWithoutReadGroupIterator( bam_fetch_iterator_t * bamIterator,
                                                                      std::string sampleName,
                                                                      utils::referenceSequencePtr_t refSequence,
                                                                      utils::timerPtr_t timer )
        : AbstractBamFileIterator( bamIterator, refSequence, timer ), m_sampleName( sampleName )
    {
    }

    BamFileWithoutReadGroupIterator::~BamFileWithoutReadGroupIterator() {}

    std::pair< std::string, readPtr_t > BamFileWithoutReadGroupIterator::getReadData()
    {
        if ( not this->hasReadData() )
        {
            throw utils::echidna_exception( "Tried to get read data without checking if there were more to get" );
        }

        readPtr_t read = std::make_shared< Read >( m_bamRecordPtr, m_refSequence );

        return std::make_pair( m_sampleName, read );
    }
}
}
