// All content Copyright (C) 2018 Genomics plc
#include <algorithm>

#include "utils/sequence.hpp"
#include "vcf/reader.hpp"
#include "caller/region.hpp"
#include "candidateVariantBank.hpp"

namespace echidna
{
namespace caller
{

    //-----------------------------------------------------------------------------------------

    namespace
    {
        struct isAlleleFrequency
        {
            bool operator()( const std::pair< std::string, std::vector< std::string > > & item )
            {
                return item.first == std::string( "AF" );
            }
        };
    }

    //-----------------------------------------------------------------------------------------

    std::vector< double > getPriorsFromInfo( const vcf::Info & info )
    {
        const auto foundPriors = std::find_if( info.begin(), info.end(), isAlleleFrequency() );

        if ( foundPriors == info.end() )
        {
            return std::vector< double >();
        }
        else
        {
            const auto & strPriors = foundPriors->second;
            std::vector< double > doublePriors( foundPriors->second.size() );

            std::transform( strPriors.begin(), strPriors.end(), doublePriors.begin(), []( const std::string & val )
                                                                                          -> double
                                                                                      {
                                                                                          return std::stod( val );
                                                                                      } );

            return doublePriors;
        }
    }

    //-----------------------------------------------------------------------------------------
}
}
