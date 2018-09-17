// All content Copyright (C) 2018 Genomics plc
#ifndef VCF_READER_HPP
#define VCF_READER_HPP

#include <utils/logging.hpp>
#include <fstream>
#include <boost/optional/optional.hpp>
#include <map>
#include <set>

#include "vcf/filterDescription.hpp"
#include "vcf/record.hpp"

#include "utils/timer.hpp"

struct shouldParseValidVCFFilterHeaderUpperCaseID;
struct shouldParseValidVCFFilterHeaderLowerCaseID;
struct shouldParseValidVCFFilterHeaderDigitsAndPunctuation;
struct shouldRaiseOnINFOHeaderType;
struct shouldRaiseOnInvalidVCFFilterHeader;
struct shouldRaiseOnFORMATHeaderType;

namespace echidna
{
namespace vcf
{
    using vcfMetaInformation_t = std::map< std::string, std::string >;

    Info parseVCFInfo( const std::string & raw_info );
}
}

#endif
