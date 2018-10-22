// All content Copyright (C) 2018 Genomics plc
#include <string>
#include <vector>

#include "caller/region.hpp"

#include <tabix/tabix.h>

namespace wecall
{
namespace io
{
    const std::string headerPrefix = "#";

    class TabixFile
    {
    public:
        explicit TabixFile( std::string filename, std::string indexFilename );

        ~TabixFile();

        std::vector< std::string > header() const { return m_headerLines; }
        std::vector< std::string > fetch( const caller::Region & region ) const;

    private:
        void readHeader();

        tabix_t * m_tabixFile;
        std::vector< std::string > m_headerLines;
    };
}
}
