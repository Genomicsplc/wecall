// All content Copyright (C) 2018 Genomics plc
#pragma once
#ifndef CALLER_PARAMS_HPP
#define CALLER_PARAMS_HPP

#include <boost/filesystem.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/optional/optional.hpp>

#include "utils/exceptions.hpp"
#include <limits>

#include "common.hpp"
#include "caller/region.hpp"
#include "vcf/field.hpp"
#include "regionUtils.hpp"

namespace echidna
{
namespace caller
{
    namespace params
    {
        namespace defaults
        {
            // -----------------------------------------------------------------------
            // Calling Params

            const bool recalibrateBaseQs = false;
            const bool outputPhasedGenotypes = true;
            const double alleleBiasThreshP = 0.9e-2;
            const double strandBiasThreshP = 1e-2;
            const double allelePlusStrandBiasThreshP = 7e-2;
            const phred_t minRMSMappingQ = 25.0;
            const double minSNPQOverDepth = 8.0;
            const double minIndelQOverDepth = 3.5;
            const phred_t minBaseQual = 25;
            const bool allowMNPCalls = false;
            const bool normalizeVariantCalls = false;
            const int minReadsToMakeCombinationClaim = 3;
            const bool turnOnLargeVariantCalls = false;

            // -----------------------------------------------------------------------
            // Private Data Params

            const std::string candidateVariantsFile;
            const std::string intermediateRecalibFileStem;

            // -----------------------------------------------------------------------
            // Private Calling Params

            const int minReadsPerVar = 2;
            const int perSamPercentReadsPerVar = 4;
            const int maxHaplotypesPerCluster = 32;
            const int minClusterDist = 20;
            const int maxClusterDist = 40;
            const int defaultBreakpointKmerSize = 7;
            const int maxBreakpointKmerSize = 70;
            const int maxClusterSize = 250;
            const int largeVariantSizeDefinition = 10;
            const int maxClusterVariantCombinations = 1024;
            const int maxClusterVariants = 30;
            const int regionPadding = maxClusterSize;  // Default to 2 * cluster distance.
            const phred_t minCallQual = 10;
            const std::string varFilterIDs = boost::algorithm::join( vcf::filter::VCFKeys, "," );
            const int badReadsWindowSize = 7;
            const phred_t minBadReadsScore = 0.0;
            const bool allVariants = false;
            const unsigned int ploidy = 2;
            const double referenceCallQualityDeltaThreshold = 0.1;

            // -----------------------------------------------------------------------
            // Filter Params

            const phred_t readMappingFilterQ = 20;
            const int baseCallFilterN = 20;
            const phred_t baseCallFilterQ = 20;
            const bool duplicatesFilter = true;
            const bool noSimilarReadsFilter = true;
            const bool noMatesFilter = true;
            const bool overlapTrim = true;
            const bool shortReadFilter = true;
            const bool shortReadTrim = true;
            const bool allowImproperPairs = true;

            // -----------------------------------------------------------------------
            // System Params

            const bool overwrite = true;

            const std::size_t maxBlockSize = 1100000;
            const std::size_t maxBlockSizeMin = 1;
            const std::size_t maxBlockSizeMax = 5000000;

            const std::size_t memLimitDefault = 1024;
            const std::size_t memLimitMin = 50;
            const std::size_t memLimitMax = 1024 * 1024;

            const std::size_t numberOfJobsDefault = 0;
            const std::size_t numberOfJobsMin = 0;
            const std::size_t numberOfJobsMax = 64;

            // -----------------------------------------------------------------------
            // Data Params

            const std::string genotypeAllelesFile = std::string();
            const std::string regions = "";
            const std::string output = "output.vcf";
            const std::string outputFormat = constants::vcf42;
            const std::string workDirDefault = "workDir";
            const bool outputRefCalls = false;
            const std::size_t maxRefCallSize = 1000;

            // -----------------------------------------------------------------------
            // Logging Params

            const int logLevel = loggingLevel::INFO;
            const std::string logFilename = constants::weCallString + ".log";
            const bool quietMode = false;
            const int verbosity = loggingLevel::INFO;
            const bool logTimings = false;

            // -----------------------------------------------------------------------
            // Private System Params

            const std::size_t biteSize = 1000;
        }

        const std::vector< std::string > allowableOutputFormats = {constants::vcf41, constants::vcf42};

        template < typename Sequence >
        std::string displayOptions( Sequence sequence )
        {
            if ( sequence.empty() )
            {
                return "";
            }
            else if ( sequence.size() == 1 )
            {
                std::stringstream ssOptions;
                ssOptions << sequence.front();
                return ssOptions.str();
            }
            else
            {
                std::stringstream ssOptions;

                for ( unsigned int i = 0; i < sequence.size() - 2; ++i )
                {
                    ssOptions << sequence[i] << ", ";
                }
                ssOptions << sequence[sequence.size() - 2] << " or " << sequence[sequence.size() - 1];

                return ssOptions.str();
            }
        }

        using boost::program_options::variables_map;
        using boost::program_options::options_description;

        class VcfOptionsLog
        {
        public:
            void addOptions( const boost::program_options::basic_parsed_options< char > & parsedOptions );
            std::string str() const;

        private:
            std::map< std::string, std::vector< std::string > > m_log;
        };

        constexpr std::size_t PARAM_HELP_DISPLAY_WIDTH = 160;

        // Template function for validating a parameter is within limits
        template < typename T >
        void validateParam( const T param,
                            const boost::optional< T > min,
                            const boost::optional< T > max,
                            const std::string & paramName )
        {
            if ( ( min.is_initialized() and param < min.get() ) or ( max.is_initialized() and param > max.get() ) )
            {
                std::stringstream msg;
                msg << "<" << paramName << "> not in acceptable range. ";
                throw utils::echidna_exception( msg.str() );
            }
        }

        // Template function for pulling in a parameter
        template < typename T >
        T getParam( const std::string & name,
                    const variables_map & optValues,
                    const boost::optional< T > min = boost::none,
                    const boost::optional< T > max = boost::none )
        {
            if ( optValues.count( name ) )
            {
                const T param = optValues[name].as< T >();
                validateParam( param, min, max, name );
                return param;
            }
            else
            {
                throw utils::echidna_exception( "Missing parameter value for " + name + "." );
            }
        }

        // Specific function for pulling in multi-value, comma-delimited string values
        std::vector< std::string > getParamList( const std::string & name, const variables_map & optValues );
        std::vector< std::string > getParamList( std::string commaSeparatedInputs );

        struct Application
        {
            Application( std::string appName,
                         std::string appVersion,
                         std::string appCommit,
                         std::string buildDate,
                         std::string appOptions )
                : m_appName( appName ),
                  m_appVersion( appVersion ),
                  m_appCommit( appCommit ),
                  m_buildDate( buildDate ),
                  m_appOptions( appOptions )
            {
            }

            std::string m_appName;
            std::string m_appVersion;
            std::string m_appCommit;
            std::string m_buildDate;
            std::string m_appOptions;
        };

        class Reduce
        {
        public:
            Reduce( const variables_map & optValues )
                : m_inputDir( getParam< std::string >( "inputDir", optValues ) ),
                  m_outputDataSink( getParam< std::string >( "output", optValues ) )
            {
            }

            Reduce( std::string inputDir, std::string outputDataSink )
                : m_inputDir( inputDir ), m_outputDataSink( outputDataSink )
            {
            }

            static options_description getOptionsDescription();

            std::string inputDir() const { return m_inputDir; }
            std::string outputDataSink() const { return m_outputDataSink; }

        private:
            std::string m_inputDir;
            std::string m_outputDataSink;
        };

        void validateReduceParams( const Reduce & reduceParams );
        void validateReduceParamsPreMap( const Reduce & reduceParams );

        class Data
        {
        public:
            Data( const variables_map & optValues, bool overwrite );

            static options_description getOptionsDescription();

            std::vector< Data > splitWorkload() const;

            std::vector< std::string > const & inputDataSources() const { return m_inputDataSources; }
            std::string outputDataSink() const { return m_outputDataSink; }
            std::string outputFormat() const { return m_outputFormat; }
            std::string workDir() const { return m_workDir; }
            std::string refFile() const { return m_refFile; }
            const partitionedRegions_t & dataRegions() const { return m_dataRegions; }
            bool outputRefCalls() const { return m_outputRefCalls; }
            std::size_t maxRefCallSize() const { return m_maxRefCallSize; }

        private:
            int64_t totalRegionLength() const;

            Data( std::vector< std::string > const & inputDataSources,
                  const std::string & outputDataSink,
                  const std::string & outputFormat,
                  const std::string & workDir,
                  const std::string & refFile,
                  const partitionedRegions_t & dataRegions,
                  const bool & outputRefCalls,
                  const std::size_t & maxRefCallSize )
                : m_inputDataSources( inputDataSources ),
                  m_outputDataSink( outputDataSink ),
                  m_outputFormat( outputFormat ),
                  m_workDir( workDir ),
                  m_refFile( refFile ),
                  m_dataRegions( dataRegions ),
                  m_outputRefCalls( outputRefCalls ),
                  m_maxRefCallSize( maxRefCallSize )
            {
            }

        private:
            std::vector< std::string > m_inputDataSources;
            std::string m_outputDataSink;
            std::string m_outputFormat;
            std::string m_workDir;
            std::string m_refFile;
            partitionedRegions_t m_dataRegions;
            bool m_outputRefCalls;
            std::size_t m_maxRefCallSize;
        };

        struct PrivateData
        {
            PrivateData( const variables_map & optValues )
                : m_candidateVariantsFile( getParam< std::string >( "candidateVariantsFile", optValues ) ),
                  m_intermediateRecalibFileStem( getParam< std::string >( "intermediateRecalibFileStem", optValues ) ),
                  m_genotypeAllelesFile( getParam< std::string >( "genotypeAllelesFile", optValues ) )
            {
                if ( this->genotypingMode() )
                {
                    //                    ECHIDNA_ERROR( ( not m_outputRefCalls ),
                    //                                   "Genotyping is incompatible with outputting reference calls."
                    //                                   );

                    ECHIDNA_ERROR( boost::filesystem::exists( this->m_genotypeAllelesFile ),
                                   "Genotype file " + this->m_genotypeAllelesFile + " does not exist" );
                    ECHIDNA_ERROR(
                        ( boost::filesystem::path( this->m_genotypeAllelesFile ).extension().string() == ".gz" ),
                        "File " + this->m_genotypeAllelesFile + " does not have .gz extension" );
                    ECHIDNA_ERROR( boost::filesystem::exists( this->m_genotypeAllelesFile + ".tbi" ),
                                   "Genotype index file " + this->m_genotypeAllelesFile + ".tbi does not exist" );
                }
            }

            static options_description getOptionsDescription();

            bool genotypingMode() const { return m_genotypeAllelesFile != defaults::genotypeAllelesFile; }
            std::string genotypeAllelesFile() const { return m_genotypeAllelesFile; }

            std::string m_candidateVariantsFile;
            std::string m_intermediateRecalibFileStem;
            std::string m_genotypeAllelesFile;
        };

        struct Logging
        {
            Logging( const variables_map & optValues )
                : m_logLevel( getParam< int >( "logLevel", optValues, 0 ) ),
                  m_logFilename( getParam< std::string >( "logFilename", optValues ) ),
                  m_quietMode( getParam< bool >( "quietMode", optValues ) ),
                  m_verbosity( getParam< int >( "verbosity", optValues ) ),
                  m_logTimings( getParam< bool >( "logTimings", optValues ) )
            {
            }

            Logging( int logLevel, std::string logFileName, bool quietMode, int verbosity, bool logTimings )
                : m_logLevel( logLevel ),
                  m_logFilename( logFileName ),
                  m_quietMode( quietMode ),
                  m_verbosity( verbosity ),
                  m_logTimings( logTimings )
            {
            }

            static options_description getOptionsDescription();

            int m_logLevel;
            std::string m_logFilename;
            bool m_quietMode;
            int m_verbosity;
            bool m_logTimings;
        };

        struct PrivateSystem
        {
            PrivateSystem( variables_map & optValues )
                : m_overwrite( getParam< bool >( "overwrite", optValues ) ),
                  m_biteSize( getParam< std::size_t >( "biteSize", optValues, 0ul ) ),
                  m_memLimit(
                      getParam< std::size_t >( "memLimit", optValues, defaults::memLimitMin, defaults::memLimitMax ) )
            {
            }

            static options_description getOptionsDescription();

            bool m_overwrite;
            std::size_t m_biteSize;
            std::size_t m_memLimit;
        };

        struct System
        {
            System( const variables_map & optValues )
                : m_maxBlockSize( getParam< std::size_t >( "maxBlockSize",
                                                           optValues,
                                                           defaults::maxBlockSizeMin,
                                                           defaults::maxBlockSizeMax ) ),
                  m_numberOfJobs( getParam< std::size_t >( "numberOfJobs",
                                                           optValues,
                                                           defaults::numberOfJobsMin,
                                                           defaults::numberOfJobsMax ) )

            {
            }

            static options_description getOptionsDescription();

            std::size_t m_maxBlockSize;
            std::size_t m_numberOfJobs;
        };

        struct Filters
        {
            Filters( const variables_map & optValues )
                : m_readMappingFilterQ( getParam< phred_t >( "readMappingFilterQ", optValues ) ),
                  m_baseCallFilterN( getParam< int >( "baseCallFilterN", optValues ) ),
                  m_baseCallFilterQ( getParam< phred_t >( "baseCallFilterQ", optValues ) ),
                  m_duplicatesFilter( getParam< bool >( "duplicatesFilter", optValues ) ),
                  m_noSimilarReadsFilter( getParam< bool >( "noSimilarReadsFilter", optValues ) ),
                  m_noMatesFilter( getParam< bool >( "noMatesFilter", optValues ) ),
                  m_overlapTrim( getParam< bool >( "overlapTrim", optValues ) ),
                  m_shortReadFilter( getParam< bool >( "shortReadFilter", optValues ) ),
                  m_shortReadTrim( getParam< bool >( "shortReadTrim", optValues ) ),
                  m_allowImproperPairs( getParam< bool >( "allowImproperPairs", optValues ) )
            {
            }

            Filters()
                : m_readMappingFilterQ( defaults::readMappingFilterQ ),
                  m_baseCallFilterN( defaults::baseCallFilterN ),
                  m_baseCallFilterQ( defaults::baseCallFilterQ ),
                  m_duplicatesFilter( defaults::duplicatesFilter ),
                  m_noSimilarReadsFilter( defaults::noSimilarReadsFilter ),
                  m_noMatesFilter( defaults::noMatesFilter ),
                  m_overlapTrim( defaults::overlapTrim ),
                  m_shortReadFilter( defaults::shortReadFilter ),
                  m_shortReadTrim( defaults::shortReadTrim ),
                  m_allowImproperPairs( defaults::allowImproperPairs )
            {
            }

            static options_description getOptionsDescription();

            phred_t m_readMappingFilterQ;
            int m_baseCallFilterN;
            phred_t m_baseCallFilterQ;
            bool m_duplicatesFilter;
            bool m_noSimilarReadsFilter;
            bool m_noMatesFilter;
            bool m_overlapTrim;
            bool m_shortReadFilter;
            bool m_shortReadTrim;
            bool m_allowImproperPairs;
        };

        struct PrivateCalling
        {
            PrivateCalling( const variables_map & optValues )
                : m_minReadsPerVar( getParam< int >( "minReadsPerVar", optValues ) ),
                  m_perSamPercentReadsPerVar( getParam< int >( "perSamPercentReadsPerVar", optValues ) ),
                  m_minBaseQual( getParam< phred_t >( "minBaseQual", optValues ) ),
                  m_maxHaplotypesPerCluster( getParam< int >( "maxHaplotypesPerCluster", optValues ) ),
                  m_minClusterDist( getParam< int >( "minClusterDist", optValues ) ),
                  m_maxClusterDist( getParam< int >( "maxClusterDist", optValues ) ),
                  m_maxClusterSize( getParam< int >( "maxClusterSize", optValues ) ),
                  m_largeVariantSizeDefinition( getParam< int >( "largeVariantSizeDefinition", optValues ) ),
                  m_maxClusterVariantCombinations( getParam< int >( "maxClusterVariantCombinations", optValues ) ),
                  m_maxClusterVariants( getParam< int >( "maxClusterVariants", optValues ) ),
                  m_regionPadding( getParam< int >( "regionPadding", optValues ) ),
                  m_varFilterIDs( getParamList( "varFilterIDs", optValues ) ),
                  m_allVariants( getParam< bool >( "allVariants", optValues ) ),
                  m_ploidy( getParam< unsigned int >( "ploidy", optValues ) ),
                  m_referenceCallQualityDeltaThreshold(
                      getParam< double >( "referenceCallQualityDeltaThreshold", optValues ) ),
                  m_normalizeVariantCalls( getParam< bool >( "normalizeVariantCalls", optValues ) ),
                  m_minReadsToMakeCombinationClaim( getParam< int >( "minReadsToMakeCombinationClaim", optValues ) ),
                  m_turnOnLargeVariantCalls( getParam< bool >( "turnOnLargeVariantCalls", optValues ) )
            {
                std::vector< std::string > unrecognisedFilterIDs = {};
                for ( const auto & varFilterID : m_varFilterIDs )
                {
                    if ( std::find( vcf::filter::VCFKeys.cbegin(), vcf::filter::VCFKeys.cend(), varFilterID ) ==
                         vcf::filter::VCFKeys.cend() )
                    {
                        unrecognisedFilterIDs.push_back( varFilterID );
                    }
                }
                ECHIDNA_ERROR( unrecognisedFilterIDs.empty(),
                               "Could not find filter ID(s): " + displayOptions( unrecognisedFilterIDs ) );
            }

            static options_description getOptionsDescription();

            int m_minReadsPerVar;
            int m_perSamPercentReadsPerVar;
            phred_t m_minBaseQual;
            int m_maxHaplotypesPerCluster;
            int m_minClusterDist;
            int m_maxClusterDist;
            int m_maxClusterSize;
            int m_largeVariantSizeDefinition;
            int m_maxClusterVariantCombinations;
            int m_maxClusterVariants;
            int m_regionPadding;
            std::vector< std::string > m_varFilterIDs;
            bool m_allVariants;
            unsigned int m_ploidy;
            double m_referenceCallQualityDeltaThreshold;
            bool m_normalizeVariantCalls;
            int m_minReadsToMakeCombinationClaim;
            bool m_turnOnLargeVariantCalls;
        };

        struct Calling
        {
            Calling( const variables_map & optValues )
                : m_allowMNPCalls( getParam< bool >( "allowMNPCalls", optValues ) ),
                  m_recalibrateBaseQs( getParam< bool >( "recalibrateBaseQs", optValues ) ),
                  m_outputPhasedGenotypes( getParam< bool >( "outputPhasedGenotypes", optValues ) ),
                  m_minAlleleBiasP( getParam< double >( "minAlleleBiasP", optValues ) ),
                  m_minStrandBiasP( getParam< double >( "minStrandBiasP", optValues ) ),
                  m_minAllelePlusStrandBiasP( getParam< double >( "minAllelePlusStrandBiasP", optValues ) ),
                  m_minRMSMappingQ( getParam< phred_t >( "minRMSMappingQ", optValues ) ),
                  m_minSNPQOverDepth( getParam< double >( "minSNPQOverDepth", optValues ) ),
                  m_minIndelQOverDepth( getParam< double >( "minIndelQOverDepth", optValues ) ),
                  m_minCallQual( getParam< phred_t >( "minCallQual", optValues ) ),
                  m_badReadsWindowSize( getParam< int >( "badReadsWindowSize", optValues, 0, 100 ) ),
                  m_minBadReadsScore( getParam< phred_t >( "minBadReadsScore", optValues, 0. ) )
            {
            }

            static options_description getOptionsDescription();

            bool m_allowMNPCalls;
            bool m_recalibrateBaseQs;
            bool m_outputPhasedGenotypes;
            double m_minAlleleBiasP;
            double m_minStrandBiasP;
            double m_minAllelePlusStrandBiasP;
            phred_t m_minRMSMappingQ;
            double m_minSNPQOverDepth;
            double m_minIndelQOverDepth;
            phred_t m_minCallQual;
            int m_badReadsWindowSize;
            phred_t m_minBadReadsScore;
        };
    }
}
}

#endif
