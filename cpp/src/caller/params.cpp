// All content Copyright (C) 2018 Genomics plc
#include <boost/filesystem/operations.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>

#include "vcf/reader.hpp"
#include "io/fastaFile.hpp"
#include "caller/params.hpp"
#include "caller/regionUtils.hpp"
#include "utils/logging.hpp"
#include "utils/exceptions.hpp"
#include "utils/partition.hpp"

namespace echidna
{
namespace caller
{
    namespace params
    {
        using boost::program_options::variables_map;
        using boost::program_options::options_description;
        using boost::program_options::value;

        void VcfOptionsLog::addOptions( const boost::program_options::basic_parsed_options< char > & parsedOptions )
        {
            for ( const auto & option : parsedOptions.options )
            {
                if ( m_log.find( option.string_key ) == m_log.end() )
                {
                    m_log[option.string_key] = std::vector< std::string >();
                    for ( const auto & value : option.value )
                    {
                        m_log[option.string_key].push_back( value );
                    }
                }
                else
                {
                    throw utils::echidna_exception( "Duplicate command line/config parameter \"--" + option.string_key +
                                                    "\". \nAborting...." );
                }
            }
        }

        std::string VcfOptionsLog::str() const
        {
            bool isFirstEntry = true;
            std::ostringstream oss;
            for ( auto option : m_log )
            {
                if ( not isFirstEntry )
                {
                    oss << "|";
                }

                oss << " " << option.first << ": ";
                for ( auto value : option.second )
                {
                    oss << value << " ";
                }
                isFirstEntry = false;
            }

            return oss.str();
        }

        std::vector< std::string > getParamList( const std::string & name, const variables_map & optValues )
        {
            if ( optValues.count( name ) )
            {
                return getParamList( optValues[name].as< std::string >() );
            }
            else
            {
                return {};
            }
        }

        std::vector< std::string > getParamList( std::string commaSeparatedInputs )
        {
            commaSeparatedInputs.erase(
                std::remove_if( commaSeparatedInputs.begin(), commaSeparatedInputs.end(), isspace ),
                commaSeparatedInputs.end() );

            std::vector< std::string > values;
            boost::split( values, commaSeparatedInputs, boost::is_any_of( "," ), boost::algorithm::token_compress_on );

            auto isEmptyString = []( const std::string & s )
            {
                return s.empty();
            };
            values.erase( std::remove_if( values.begin(), values.end(), isEmptyString ), values.end() );

            return values;
        }

        void validateReduceParams( const Reduce & reduceParams )
        {
            ECHIDNA_LOG( INFO, "Validating and reducing params" );
            ECHIDNA_ERROR( boost::filesystem::exists( reduceParams.inputDir() ),
                           "Working dir: " + reduceParams.inputDir() + " does not exist" );
            ECHIDNA_ERROR( boost::filesystem::is_directory( reduceParams.inputDir() ),
                           "Working dir: " + reduceParams.inputDir() + " is not a directory" );
        }
        void validateReduceParamsPreMap( const Reduce & reduceParams )
        {
            if ( boost::filesystem::exists( reduceParams.inputDir() ) )
            {
                ECHIDNA_ERROR( boost::filesystem::is_empty( reduceParams.inputDir() ),
                               "Working dir: " + reduceParams.inputDir() + " is not empty" );
            }
        }

        Data::Data( const variables_map & optValues, bool overwrite )
            : m_inputDataSources( getParamList( "inputs", optValues ) ),
              m_outputDataSink( getParam< std::string >( "output", optValues ) ),
              m_outputFormat( getParam< std::string >( "outputFormat", optValues ) ),
              m_workDir( getParam< std::string >( "workDir", optValues ) ),
              m_refFile( getParam< std::string >( "refFile", optValues ) ),
              m_outputRefCalls( getParam< bool >( "outputRefCalls", optValues ) ),
              m_maxRefCallSize( getParam< std::size_t >( "maxRefCallSize", optValues ) )
        {

            ECHIDNA_ERROR( ( std::find( allowableOutputFormats.cbegin(), allowableOutputFormats.cend(),
                                        m_outputFormat ) != allowableOutputFormats.cend() ),
                           std::string( "output file format must be " + displayOptions( allowableOutputFormats ) ) );

            for ( auto const & inputDataSource : m_inputDataSources )
            {
                ECHIDNA_ERROR( boost::filesystem::exists( inputDataSource ),
                               "File " + inputDataSource + " does not exist" );
                ECHIDNA_ERROR( ( boost::filesystem::path( inputDataSource ).extension().string() == ".bam" ),
                               "File " + inputDataSource + " does not have .bam extension" );
                ECHIDNA_ERROR( boost::filesystem::exists( inputDataSource + ".bai" ),
                               "Index file " + inputDataSource + ".bai does not exist" );
            }
            if ( not overwrite )
            {
                ECHIDNA_ERROR( not boost::filesystem::exists( m_outputDataSink ),
                               m_outputDataSink + " already exists" );
            }
            ECHIDNA_ERROR( boost::filesystem::exists( m_refFile ), "File " + m_refFile + " does not exist" );
            ECHIDNA_ERROR( ( boost::filesystem::path( m_refFile ).extension().string() == ".fa" ),
                           "File " + m_refFile + " does not have .fa extension" );

            const auto expectedFastaIndexFile = io::fastaIndexFileName( m_refFile );
            ECHIDNA_ERROR( boost::filesystem::exists( expectedFastaIndexFile ),
                           "Index file " + expectedFastaIndexFile + " does not exist" );

            m_dataRegions = DataRegionsBuilder( getParamList( "regions", optValues ),
                                                io::FastaIndex( expectedFastaIndexFile ) ).build();
        }

        options_description Reduce::getOptionsDescription()
        {
            options_description options( "Reduce Parameters", PARAM_HELP_DISPLAY_WIDTH );

            options.add_options()( "inputDir", value< std::string >()->required(),
                                   "input directory containing input VCF files" )(
                "output", value< std::string >()->required(), "output VCF data file name" );
            return options;
        }

        // clang-format off
        options_description Data::getOptionsDescription()
        {
            options_description options("Data Parameters", PARAM_HELP_DISPLAY_WIDTH);

            options.add_options()
                ("inputs", value<std::string>()->required(), "comma separated list of input BAM data file names")
                ("refFile", value<std::string>()->required(), "reference genome file")
                ("regions", value<std::string>()->default_value(defaults::regions), "regions to process -- comma separated list of bed files or of chroms or chrom:start-end's.")
                ("output", value<std::string>()->default_value(defaults::output), "output file name")
                ("outputFormat", value<std::string>()->default_value(defaults::outputFormat), std::string("output file format (" + displayOptions(allowableOutputFormats) + ")").c_str())
                ("workDir", value<std::string>()->default_value(defaults::workDirDefault), "intermediate files directory (for parallel runs only)")
                ("outputRefCalls", value<bool>()->default_value(defaults::outputRefCalls)->implicit_value(true), "if specified, output reference as well as variant calls")
                ("maxRefCallSize", value<std::size_t>()->default_value(defaults::maxRefCallSize), "maximum size of individual reference calls")
                ;

            return options;
        }
        // clang-format on

        namespace
        {
            void validateAndCreateWorkingDir( const std::string & dir )
            {
                if ( boost::filesystem::exists( dir ) )
                {
                    ECHIDNA_ERROR( boost::filesystem::is_directory( dir ),
                                   "Working dir: " + dir + " is not a directory" );
                }
                else
                {
                    bool success = boost::filesystem::create_directory( dir );
                    ECHIDNA_ERROR( success, "Could not create work directory: " + dir );
                }
            }
        }

        std::vector< Data > Data::splitWorkload() const
        {
            std::vector< Data > vecData;

            validateAndCreateWorkingDir( m_workDir );

            boost::format intermediateFileNameFormat( "%05d.vcf" );
            ECHIDNA_ERROR( ( m_dataRegions.size() < 99999 ),
                           constants::weCallString + " called with too many regions. Max=99999" );

            for ( std::size_t i = 0; i < m_dataRegions.size(); ++i )
            {
                boost::filesystem::path outputDataSink = boost::filesystem::path( m_workDir );
                outputDataSink /= ( intermediateFileNameFormat % i ).str();
                ECHIDNA_LOG( DEBUG, outputDataSink.string() );

                ECHIDNA_ERROR( ( not boost::filesystem::exists( outputDataSink ) ),
                               "output data sink " + outputDataSink.string() + " already exist" );

                vecData.push_back( Data( m_inputDataSources, outputDataSink.string(), m_outputFormat, m_workDir,
                                         m_refFile, {m_dataRegions[i]}, m_outputRefCalls, m_maxRefCallSize ) );
            }

            std::sort( vecData.begin(), vecData.end(), []( const Data & first, const Data & second )
                       {
                           return first.totalRegionLength() > second.totalRegionLength();
                       } );

            return vecData;
        }

        int64_t Data::totalRegionLength() const
        {
            int64_t totalLength = 0;
            for ( const auto & vecRegions : m_dataRegions )
            {
                for ( const auto & region : vecRegions )
                {
                    totalLength += region.size();
                }
            }

            return totalLength;
        }

        std::string logLevelsString()
        {
            std::stringstream message;
            message << "(";
            for ( std::size_t i = 0; i < numberLogLevels; ++i )
            {
                message << i << ": " << logLabels[i];
                if ( i != numberLogLevels - 1 )
                {
                    message << ", ";
                }
            }
            message << ")";
            return message.str();
        }

        // clang-format off
        options_description Logging::getOptionsDescription()
        {
            options_description options("Logging parameters", PARAM_HELP_DISPLAY_WIDTH);

            options.add_options()
                ("logLevel", value<int>()->default_value(defaults::logLevel), ("logging level for log file " + logLevelsString()).c_str())
                ("logFilename", value<std::string>()->default_value(defaults::logFilename), "path to file where logging messages will be written")
                ("quietMode", value<bool>()->default_value(defaults::quietMode)->implicit_value(true), "no stream output except errors")
                ("verbosity", value<int>()->default_value(defaults::verbosity), ("verbosity level for stdout " + logLevelsString()).c_str())
                ("logTimings", value<bool>()->default_value(defaults::logTimings), "output timings log messages")
                ;

            return options;
        }

        options_description PrivateSystem::getOptionsDescription()
        {
            options_description options("Config only parameters", PARAM_HELP_DISPLAY_WIDTH);

            options.add_options()
                ("overwrite", value<bool>()->default_value(defaults::overwrite), "allow overwrite of output vcf")
                ("biteSize", value <std::size_t>()->default_value(defaults::biteSize), "size of read data buffer increment (in bases)")
                    ("memLimit", value <std::size_t>()->default_value(defaults::memLimitDefault))
                ;

            return options;
        }

        options_description System::getOptionsDescription()
        {
            options_description options("System parameters", PARAM_HELP_DISPLAY_WIDTH);

            const std::string block_message = "maximum block size -- larger regions will be broken up into blocks. Must be within the range " +
                std::to_string(defaults::maxBlockSizeMin) + "-" + std::to_string(defaults::maxBlockSizeMax) + " inclusive.";

            const std::string jobs_message = "number of jobs to run simultaneously. Must be within the range " +
                std::to_string(defaults::numberOfJobsMin) + "-" + std::to_string(defaults::numberOfJobsMax) +
                " inclusive.";

            options.add_options()
                ("maxBlockSize", value <std::size_t>()->default_value(defaults::maxBlockSize), block_message.c_str())
                ("numberOfJobs", value <std::size_t>()->default_value(defaults::numberOfJobsDefault), jobs_message.c_str())
                ;

            return options;
        }

        options_description Filters::getOptionsDescription()
        {
            options_description options("Read filter parameters", PARAM_HELP_DISPLAY_WIDTH);

            options.add_options()
                ("readMappingFilterQ", value<phred_t>()->default_value(defaults::readMappingFilterQ), "read mapping filter - minimum quality (as a phred score)")
                ("baseCallFilterN", value<int>()->default_value(defaults::baseCallFilterN), "base call filter - minimum number of 'good' base calls per read")
                ("baseCallFilterQ", value<phred_t>()->default_value(defaults::baseCallFilterQ), "base call filter - minimum quality for a base call to be 'good' (as a phred score)")
                ("duplicatesFilter", value<bool>()->default_value(defaults::duplicatesFilter)->implicit_value(true), "filter out duplicate reads")
                ("noSimilarReadsFilter", value<bool>()->default_value(defaults::noSimilarReadsFilter)->implicit_value(true), "filter out reads with matching position, sequence, cigar & quality")
                ("noMatesFilter", value<bool>()->default_value(defaults::noMatesFilter)->implicit_value(true), "filter out reads whose mate is unmapped")
                ("overlapTrim", value<bool>()->default_value(defaults::overlapTrim)->implicit_value(true), "where fragment is < 2 read lengths, trim the area of overlap on one read by setting base quality to 0 (avoids double counting)")
                ("shortReadFilter", value<bool>()->default_value(defaults::shortReadFilter)->implicit_value(true), "filter out reads where fragment < 1 read length")
                ("shortReadTrim", value<bool>()->default_value(defaults::shortReadTrim)->implicit_value(true), "trim the adapter sequence from reads where fragment < 1 read length (only relevant if shortReadFilter is off)")
                ("allowImproperPairs", value<bool>()->default_value(defaults::allowImproperPairs)->implicit_value(true), "If true, filter out read-pairs which are not labelled as proper pairs (in the BAM record bit flag)")
                ;

            return options;
        }

        options_description PrivateCalling::getOptionsDescription()
        {
            options_description options("Private calling parameters", PARAM_HELP_DISPLAY_WIDTH);
            options.add_options()
                ("minReadsPerVar", value<int>()->default_value(defaults::minReadsPerVar), "Minimum number of reads needed to support a call")
                ("perSamPercentReadsPerVar", value<int>()->default_value(defaults::perSamPercentReadsPerVar), "Must have sample with percentage read coverage greater than this to make a call")
                ("minBaseQual", value<phred_t>()->default_value(defaults::minBaseQual), "Minimum base quality for base mismatches to be considered as potential SNPs")
                ("maxHaplotypesPerCluster", value<int>()->default_value(defaults::maxHaplotypesPerCluster), "Maximum number of haplotypes considered in one cluster")
                ("minClusterDist", value<int>()->default_value(defaults::minClusterDist), "Variants closer together than this must belong to the same cluster regardless of cluster size.")
                ("maxClusterDist", value<int>()->default_value(defaults::maxClusterDist), "Variants closer together than this must belong to the same cluster only if cluster size is small or cluster variants small.")
                ("maxClusterSize", value<int>()->default_value(defaults::maxClusterSize), "Upper limit on the size (in the reference) of a cluster")
                ("largeVariantSizeDefinition", value<int>()->default_value(defaults::largeVariantSizeDefinition), "Variants with ref/alt longer than this are considered to be large.")
                ("largeVariantClusterThreshold", value<int>()->default_value(defaults::largeVariantClusterThreshold), "A variant cluster is considered as a large variant cluster when it contains a variant which is longer than this threshold.")
                ("maxClusterVariantCombinations", value<int>()->default_value(defaults::maxClusterVariantCombinations), "Upper limit on the number of variants combinations merged into same cluster.")
                ("maxClusterVariants", value<int>()->default_value(defaults::maxClusterVariants), "Upper limit on the number of variants merged into same cluster.")
                ("regionPadding", value<int>()->default_value(defaults::regionPadding), "Internally pad regions by this amount whilst modelling.")
                ("varFilterIDs", value<std::string>()->default_value(defaults::varFilterIDs), "Comma separated list of variant filter IDs")
                ("allVariants", value<bool>()->default_value(defaults::allVariants), "Flag to control whether uncalled variants are output")
                ("ploidy", value<unsigned int>()->default_value(defaults::ploidy), "Ploidy used in variant calling model")
                ("referenceCallQualityDeltaThreshold", value<double>()->default_value(defaults::referenceCallQualityDeltaThreshold), "Threshold of phred-scaled quality change to start new reference call block.")
                ("normalizeVariantCalls", value<bool>()->default_value(defaults::normalizeVariantCalls), "Enable/Disable Variant representation normalization")
                ("minReadsToMakeCombinationClaim", value<int>()->default_value(defaults::minReadsToMakeCombinationClaim), "Minimum number of reads to make combination claim")
                ;

            return options;
        }

        options_description PrivateData::getOptionsDescription()
        {
            options_description options("Private calling parameters", PARAM_HELP_DISPLAY_WIDTH);
            options.add_options()
                ("candidateVariantsFile", value<std::string>()->default_value(defaults::candidateVariantsFile.c_str()), "ad-hoc candidate variants file")
                ("intermediateRecalibFileStem", value<std::string>()->default_value(defaults::intermediateRecalibFileStem.c_str()), "Intermediate Sam File stem of reads with recalibrated base qualities. Outputs stem_{sample_name}.sam")
                ("genotypeAllelesFile", value<std::string>()->default_value(defaults::genotypeAllelesFile), "variants to be genotyped")
                ;

            return options;
        }

        options_description Calling::getOptionsDescription()
        {
            options_description options("Calling parameters", PARAM_HELP_DISPLAY_WIDTH);
            options.add_options()
                ("allowMNPCalls", value<bool>()->default_value(defaults::allowMNPCalls), "enable/disable MNP calling")
                ("recalibrateBaseQs", value<bool>()->default_value(defaults::recalibrateBaseQs)->implicit_value(true), "if specified, run error recalibration on read base qualities")
                ("outputPhasedGenotypes", value<bool>()->default_value(defaults::outputPhasedGenotypes), "output genotype phase information.")
                ("minAlleleBiasP", value<double>()->default_value( defaults::alleleBiasThreshP, utils::toString( defaults::alleleBiasThreshP ) ), ("minimum allele bias P-value (INFO::" + vcf::info::ABPV_key + "), below which variants will be filtered").c_str())
                ("minStrandBiasP", value<double>()->default_value( defaults::strandBiasThreshP, utils::toString( defaults::strandBiasThreshP ) ), ("minimum strand bias P-value (INFO::" + vcf::info::SBPV_key + "), below which variants will be filtered").c_str())
                ("minAllelePlusStrandBiasP", value<double>()->default_value( defaults::allelePlusStrandBiasThreshP, utils::toString( defaults::allelePlusStrandBiasThreshP ) ), ("minimum allele plus strand bias P-value (INFO::" + vcf::info::SBPV_key + " + INFO::" + vcf::info::ABPV_key + "), below which variants will be filtered").c_str())
                ("minRMSMappingQ", value<phred_t>()->default_value(defaults::minRMSMappingQ, utils::toString( defaults::minRMSMappingQ )  ), ("minimum phred-scaled root mean square mapping quality of supporting reads (INFO::" + vcf::info::MQ_key + "), below which variants will be filtered").c_str())
                ("minSNPQOverDepth", value<double>()->default_value(defaults::minSNPQOverDepth), ("minimum phred-scaled variant quality over depth (INFO::" + vcf::info::QD_key + ") for SNPs or MNPs, below which variants will be filtered").c_str())
                ("minIndelQOverDepth", value<double>()->default_value(defaults::minIndelQOverDepth), ("minimum phred-scaled variant quality over depth (INFO::" + vcf::info::QD_key + ") for Indels, below which variants will be filtered").c_str())
                ("minCallQual", value<phred_t>()->default_value(defaults::minCallQual), ("minimum phred-scaled quality for (INFO::" + vcf::info::PP_key + "), below which variants will be filtered").c_str())
                ("minBadReadsScore", value<phred_t >()->default_value(defaults::minBadReadsScore), ("minimum median of the min base qualities in bad reads window (INFO::" + vcf::info::BR_key + "), below which variants will be filtered").c_str())
                ("badReadsWindowSize", value<int>()->default_value(defaults::badReadsWindowSize), "window size around variant in which poor base quality is considered")
                ("turnOnLargeVariantCalls", value<bool>()->default_value(defaults::turnOnLargeVariantCalls), "Set to turn on Large variant calling")
                ;

            return options;
        }
        // clang-format on
    }
}
}
