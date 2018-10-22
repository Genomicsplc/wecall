// All content Copyright (C) 2018 Genomics plc
#include "readrecalibration/intermediateOutputWriter.hpp"

#include <vector>
#include <boost/filesystem/operations.hpp>
#include <fstream>

#include "io/bamFile.hpp"
#include "utils/exceptions.hpp"

namespace wecall
{
namespace corrector
{

    IntermediateOutputWriter::IntermediateOutputWriter( const std::vector< std::string > & inputBams,
                                                        std::string outputFileStem )
        : m_outputFileStem( outputFileStem ), m_writeOutputFile( not m_outputFileStem.empty() )
    {
        if ( m_writeOutputFile )
        {
            for ( const auto & inputBamFilename : inputBams )
            {
                auto inputBam = io::BamFile( inputBamFilename );
                auto sampleNames = inputBam.getSampleNames();
                if ( sampleNames.size() != 1 )
                {
                    throw utils::wecall_exception( "Only one sample per bam file is currently supported." );
                }
                auto outputFilename = this->sampleFilename( sampleNames.front() );

                m_sampleNameToFileMap[sampleNames.front()] = outputFilename;

                this->writeSamHeader( outputFilename, inputBam.bamHeader().text );
            }
        }
    }

    //-----------------------------------------------------------------------------------------

    void IntermediateOutputWriter::writeSamHeader( std::string outputFilename, std::string headerText ) const
    {
        std::ofstream outputSam;
        outputSam.open( outputFilename, std::ofstream::out | std::ofstream::trunc );
        outputSam << headerText;
        outputSam.close();
    }

    //-----------------------------------------------------------------------------------------

    std::string addExclamationMark( std::string strWithoutExclamationMark )
    {
        std::string result = strWithoutExclamationMark;
        std::transform( result.begin(), result.end(), result.begin(),
                        []( const char & charWithoutExclamationMark ) -> char
                        {
                            return static_cast< char >( charWithoutExclamationMark + 33 );
                        } );

        return result;
    }

    //-----------------------------------------------------------------------------------------

    void IntermediateOutputWriter::writeReads( const io::perSampleRegionsReads_t & readRangesPerSample,
                                               std::string contig ) const
    {
        if ( m_writeOutputFile )
        {
            for ( const auto & sampleNameReadRange : readRangesPerSample )
            {
                auto sampleName = sampleNameReadRange.first;
                auto filename = m_sampleNameToFileMap.at( sampleName );

                std::ofstream outputSam;
                outputSam.open( filename, std::ofstream::out | std::ofstream::app );

                for ( const auto & read : sampleNameReadRange.second )
                {
                    // WECALL_LOG(SUPER_DEBUG, "Outputting read with quality: " << read->getQualities());

                    outputSam << read.getQName() << "\t";
                    outputSam << read.getFlag() << "\t";
                    outputSam << contig << "\t";
                    outputSam << ( 1 + read.getStartPos() ) << "\t";
                    outputSam << read.getMappingQuality() << "\t";
                    outputSam << read.cigar().toString() << "\t";
                    outputSam << "="
                              << "\t";
                    outputSam << ( 1 + read.getMateStartPos() ) << "\t";
                    outputSam << read.getInsertSize() << "\t";
                    outputSam << read.sequence() << "\t";
                    outputSam << addExclamationMark( read.getQualities() ) << std::endl;
                }
                outputSam.close();
            }
        }
    }

    //-----------------------------------------------------------------------------------------

    std::string IntermediateOutputWriter::sampleFilename( std::string sampleName ) const
    {
        return m_outputFileStem + "_" + sampleName + ".sam";
    }

}  // namespace corrector
}  // namespace wecall
