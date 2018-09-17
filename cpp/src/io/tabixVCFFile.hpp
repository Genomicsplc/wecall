// All content Copyright (C) 2018 Genomics plc
#include <string>
#include <vector>
#include "io/tabixFile.hpp"
#include "vcf/record.hpp"
#include "vcf/filterDescription.hpp"
#include "caller/region.hpp"

namespace echidna
{
namespace io
{
    using vcfMetaInformation_t = std::map< std::string, std::string >;
    class TabixVCFFile
    {
    public:
        explicit TabixVCFFile( std::string filename, std::string indexFilename );

        ~TabixVCFFile();

        std::vector< vcf::Record > fetch( const caller::Region & region ) const;

        static bool containsFilterId( const std::set< vcf::FilterDesc > & filterDescs, const std::string & filterId );
        static std::pair< std::string, std::string > parseMetaInfoLine( std::string line );
        static vcf::FilterDesc parseFilterHeaderLine( const std::string & line );

    private:
        void readHeader( std::string filename );

        TabixFile m_tabixFile;

        vcfMetaInformation_t m_metaInformation;
        std::set< vcf::FilterDesc > m_filterDescs;
    };
}
}
