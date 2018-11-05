// All content Copyright (C) 2018 Genomics plc
#ifndef BAM_FILE_ITERATOR_HPP
#define BAM_FILE_ITERATOR_HPP

#include "utils/timer.hpp"
#include "common.hpp"
#include "io/read.hpp"

namespace wecall
{
namespace io
{

    class AbstractBamFileIterator
    {
    public:
        AbstractBamFileIterator( bam_fetch_iterator_t * bamIterator,
                                 utils::referenceSequencePtr_t refSequence,
                                 utils::timerPtr_t timer );
        virtual ~AbstractBamFileIterator();

        void next();
        bool hasReadData();
        virtual std::pair< std::string, readPtr_t > getReadData() = 0;

        utils::ScopedTimerTrigger timer() { return utils::ScopedTimerTrigger( m_timer ); }

    protected:
        bam1_t * m_bamRecordPtr;
        utils::referenceSequencePtr_t m_refSequence;

    private:
        bam_fetch_iterator_t * m_bamIterator;
        utils::timerPtr_t m_timer;
    };

    class BamFileIterator : public AbstractBamFileIterator
    {
    public:
        BamFileIterator( bam_fetch_iterator_t * bamIterator,
                         std::map< std::string, std::string > samplesByID,
                         utils::referenceSequencePtr_t refSequence,
                         utils::timerPtr_t timer );
        ~BamFileIterator();

        std::pair< std::string, readPtr_t > getReadData() override;

    private:
        const std::map< std::string, std::string > m_samplesByID;
    };

    class BamFileWithoutReadGroupIterator : public AbstractBamFileIterator
    {
    public:
        BamFileWithoutReadGroupIterator( bam_fetch_iterator_t * bamIterator,
                                         std::string sampleName,
                                         utils::referenceSequencePtr_t refSequence,
                                         utils::timerPtr_t timer );
        ~BamFileWithoutReadGroupIterator();

        std::pair< std::string, readPtr_t > getReadData() override;

    private:
        const std::string m_sampleName;
    };

    using bamFileIteratorPtr_t = std::shared_ptr< AbstractBamFileIterator >;
}
}

#endif
