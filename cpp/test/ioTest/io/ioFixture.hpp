// All content Copyright (C) 2018 Genomics plc
#ifndef IO_FIXTURE_HPP
#define IO_FIXTURE_HPP

#include "io/bamFile.hpp"
#include "io/bedFile.hpp"
#include "io/fastaFile.hpp"

#include "ioTest/utils/environment.hpp"

#include <vector>
#include <string>

namespace echidna
{
namespace test
{
    namespace fs = boost::filesystem;

    inline void writeFile( const std::string & filename, const std::string & content )
    {
        std::ofstream myfile( filename );
        if ( myfile.is_open() )
        {
            myfile << content;
            myfile.close();
        }
    }

    struct BedFileFixture
    {

        BedFileFixture()
        {
            tempDir = fs::temp_directory_path();

            bedNames.push_back( tempDir.generic_string() + "/00.bed" );
            bedNames.push_back( tempDir.generic_string() + "/01.bed" );
            bedNames.push_back( tempDir.generic_string() + "/02.bed" );
            bedNames.push_back( tempDir.generic_string() + "/03.bed" );

            writeFile( bedNames[0], "chr1\t1\t2\n" );
            writeFile( bedNames[1],
                       "chr1\t2\t3\n"
                       "chr2\t5\t7\n" );
            writeFile( bedNames[2],
                       "chr7\t0\t10\tPos1\t0\t+\t127471196\t127472363\t255,0,0\n"
                       "chr7\t20\t25\tPos2\t0\t+\t127472363\t127473530\t255,0,0\n"
                       "chr7\t30\t50\tPos3\t0\t+\t127473530\t127474697\t255,0,0\n" );
            writeFile( bedNames[3],
                       "browser position chr7:127471196-127495720\n"
                       "browser hide all\n"
                       "track name=\"ColorByStrandDemo\" description=\"Color by strand demonstration\" visibility=2 "
                       "colorByStrand=\"255,0,0 0,0,255\"\n"
                       "chr7\t0\t10\tPos1\t0\t+\t127471196\t127472363\t255,0,0\n"
                       "chr7\t20\t25\tPos2\t0\t+\t127472363\t127473530\t255,0,0\n"
                       "gi|734691289|gb|KN707645.1|\t30\t50\tPos3\t0\t+\t127473530\t127474697\t255,0,0\n" );
            std::sort( bedNames.begin(), bedNames.end() );

            for ( const auto & filename : bedNames )
            {
                bedFiles.emplace_back( new echidna::io::BedFile( filename ) );
            }
        }

        ~BedFileFixture()
        {
            for ( const auto & filename : bedNames )
            {
                fs::remove( filename.c_str() );
            }
        }

        fs::path tempDir;
        std::vector< std::string > bedNames;
        std::vector< std::unique_ptr< echidna::io::BedFile > > bedFiles;
    };

    struct FastaFileFixture
    {
        FastaFileFixture()
            : refFilename( fs::temp_directory_path().generic_string() + "/test.fa" ),
              indexFilename( fs::temp_directory_path().generic_string() + "/test.fa.fai" )
        {
            writeFile( indexFilename,
                       "1\t240\t46\t60\t61\n"
                       "2\t210\t336\t60\t61\n" );

            writeFile( refFilename,
                       ">1 dna:chromosome chromosome:GRCh37:1:1:240:1\n"
                       "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"
                       "CGCAGGCGCAGAGACACATGCTACCGCGTCCAGGGGTGGAGGCGTGGCGCAGGCGCAGAG\n"
                       "AGGCGCACCGCGCCGGCGCAGGCGCAGAGACACATGCTAGCGCGTCCAGGGGTGGAGGCG\n"
                       "TGGCGCAGGCGCAGAGACGCAAGCCTACGGGCGGGGGTTGGGGGGGCGTGTGTTGCAGGA\n"
                       ">2 dna:chromosome chromosome:GRCh37:2:1:210:1\n"
                       "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"
                       "CCACAATGCATTTGTCAAAATATGCAGAATTTTACAGCCATATGGTTAGAGCAAACTCTA\n"
                       "TTCAAATTAAATAAAATTACTCAGGATGTGGAGTATCCCAGGACAGAATACATCATGTGA\n"
                       "AAAAGCATTTATGCTACAAATTACTATGGA\n" );
            refFiles.emplace_back( new echidna::io::FastaFile( refFilename ) );
        }

        ~FastaFileFixture() { fs::remove( refFilename.c_str() ); }

        std::string refFilename;
        std::string indexFilename;
        std::vector< std::unique_ptr< echidna::io::FastaFile > > refFiles;
    };

    struct FastaIndexFileFixture
    {
        FastaIndexFileFixture() : indexFilename( fs::temp_directory_path().generic_string() + "/v01.fa.fai" )
        {
            writeFile( indexFilename,
                       "1\t249250621\t52\t60\t61\n"
                       "2\t243199373\t253404903\t60\t61\n"
                       "3\t198022430\t500657651\t60\t61\n"
                       "4\t191154276\t701980507\t60\t61\n"
                       "5\t180915260\t896320740\t60\t61\n"
                       "6\t171115067\t1080251307\t60\t61\n"
                       "7\t159138663\t1254218344\t60\t61\n"
                       "8\t146364022\t1416009371\t60\t61\n"
                       "9\t141213431\t1564812846\t60\t61\n"
                       "10\t135534747\t1708379889\t60\t61\n"
                       "11\t135006516\t1846173603\t60\t61\n"
                       "12\t133851895\t1983430282\t60\t61\n"
                       "13\t115169878\t2119513096\t60\t61\n"
                       "14\t107349540\t2236602526\t60\t61\n"
                       "15\t102531392\t2345741279\t60\t61\n"
                       "16\t90354753\t2449981581\t60\t61\n"
                       "17\t81195210\t2541842300\t60\t61\n"
                       "18\t78077248\t2624390817\t60\t61\n"
                       "19\t59128983\t2703769406\t60\t61\n"
                       "20\t63025520\t2763883926\t60\t61\n"
                       "21\t48129895\t2827959925\t60\t61\n"
                       "22\t51304566\t2876892038\t60\t61\n"
                       "X\t155270560\t2929051733\t60\t61\n"
                       "Y\t59373566\t3086910193\t60\t61\n"
                       "MT\t16569\t3147273397\t70\t71\n"
                       "GL000207.1\t4262\t3147290264\t60\t61\n"
                       "GL000192.1\t547496\t3152949897\t60\t61\n"
                       "NC_007605\t171823\t3153506529\t60\t61\n"
                       "hs37d5\t35477943\t3153681224\t60\t61\n" );
            fastaIndices.emplace_back( new echidna::io::FastaIndex( indexFilename ) );
        }

        ~FastaIndexFileFixture() { fs::remove( indexFilename.c_str() ); }

        std::string indexFilename;
        std::vector< std::unique_ptr< echidna::io::FastaIndex > > fastaIndices;
    };
}
}

#endif
