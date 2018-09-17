// All content Copyright (C) 2018 Genomics plc
#ifndef INTERMEDIATE_OUTPUT_WRITER_HPP
#define INTERMEDIATE_OUTPUT_WRITER_HPP

#include <vector>
#include <string>

#include "io/readRange.hpp"

namespace echidna
{
namespace corrector
{
    class IntermediateOutputWriter
    {

    public:
        IntermediateOutputWriter( const std::vector< std::string > & inputBams, std::string outputFileStem );

        void writeReads( const io::perSampleRegionsReads_t & readRangesPerSample, std::string contig ) const;

    private:
        std::string sampleFilename( std::string sampleName ) const;
        void writeSamHeader( std::string outputFilename, std::string headerText ) const;

    private:
        const std::string m_outputFileStem;
        std::map< std::string, std::string > m_sampleNameToFileMap;
        const bool m_writeOutputFile;
    };

}  // namespace corrector
}  // namespace echidna

#endif  // INTERMEDIATE_OUTPUT_WRITER_HPP
