// All content Copyright (C) 2018 Genomics plc
#ifndef BEDFILE_HPP
#define BEDFILE_HPP

#include <string>
#include <vector>

#include "utils/timer.hpp"
#include "caller/region.hpp"
#include "io/fastaFile.hpp"

namespace echidna
{
namespace io
{
    class BedFile
    {
    public:
        explicit BedFile( std::string fileName );

        caller::regions_t getRegions() const;

    private:
        std::vector< std::string > parseRegionLine( std::string line, bool insideHeader, size_t lineNumber ) const;

        caller::regions_t readRegionsFromStream( std::istream & inStream ) const;

    private:
        std::string m_filename;
        utils::timerPtr_t m_timer;
    };
}
}

#endif
