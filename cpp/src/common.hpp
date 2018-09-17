// All content Copyright (C) 2018 Genomics plc
#ifndef GENOMICSLTD_COMMON_HPP
#define GENOMICSLTD_COMMON_HPP

// Commonly used types with common implementations across the baseline

#include <iostream>
#include <limits>
#include <cmath>
#include "utils/exceptions.hpp"

// useful macro for printing location in source

#define STRING( arg ) STRING2( arg )
#define STRING2( arg ) #arg
#define WHERE() __FILE__ ":" STRING( __LINE__ )

// Genomics specific statistical measure of quality - the "PHRED score"

using phred_t = double;

// Logging severity levels
enum loggingLevel
{
    FATAL = 0,
    ERROR = 1,
    WARNING = 2,
    INFO = 3,
    TIMING = 4,
    DEBUG = 5,
    SUPER_DEBUG = 6
};

static const char * logLabels[] = {"FATAL", "ERROR", "WARNING", "INFO", "TIMING", "DEBUG", "SUPER_DEBUG"};

constexpr std::size_t numberLogLevels = sizeof( logLabels ) / sizeof( *logLabels );
// Function to convert log levels to strings
inline std::ostream & operator<<( std::ostream & strm, loggingLevel level )
{

    if ( static_cast< std::size_t >( level ) < numberLogLevels )
    {
        strm << logLabels[level];
    }
    else
    {
        strm << static_cast< int >( level );
    }

    return strm;
}

inline std::size_t int64_to_sizet( int64_t v )
{
    if ( v < 0 )
    {
        throw echidna::utils::echidna_exception( "Cannot convert int64_t " + std::to_string( v ) + " to size_t." );
    }
    return static_cast< std::size_t >( v );
}

inline std::size_t int32_to_sizet( int32_t v )
{
    if ( v < 0 )
    {
        throw echidna::utils::echidna_exception( "Cannot convert int32_t " + std::to_string( v ) + " to size_t." );
    }
    return static_cast< std::size_t >( v );
}

inline std::size_t long_to_sizet( long v )
{
    if ( v < 0 )
    {
        throw echidna::utils::echidna_exception( "Cannot convert long " + std::to_string( v ) + " to size_t." );
    }
    return static_cast< std::size_t >( v );
}

inline std::string substr( const std::string & input, int64_t pos, int64_t length )
{
    return input.substr( int64_to_sizet( pos ), int64_to_sizet( length ) );
}

// General (not class specific) constants that may be used across the baseline

namespace constants
{
constexpr double likelihoodCap = -300.0;  /// Arbitrary cap close to min float value
constexpr double maxPhredScore = 3000.0;
constexpr char minAllowedQualityScore = 0x02;  // The alignment code doesn't like q-scores below this
constexpr double minVariantPrior = 10e-15;
constexpr char gapChar = 'N';
const std::string vcfUnphasedGenotypeDeliminator = "/";
const std::string vcfPhasedGenotypeDeliminator = "|";
const std::string vcf41 = "VCF4.1";
const std::string vcf42 = "VCF4.2";
const std::string weCallString = "weCall";
constexpr int needlemanWunschPadding = 8;
constexpr int bamFetchRegionPadding = 100;
constexpr double strandBiasFilterBeta = 20.0;
constexpr double alleleBiasFilterAlpha = 20.0;
constexpr double alleleBiasFilterBeta = 20.0;
constexpr short gapExtendPenalty = 3;
constexpr short nucleotidePrior = 0;
const std::string vcfUnknownValue = ".";
constexpr double unknownValue = std::numeric_limits< double >::quiet_NaN();
const std::string vcfRefAltValue = "<NON_REF>";
const auto vcfRecordColumnSeparator = "\t";

constexpr double thresholdForReadSupportingVariant = 0.5;
}

#endif
