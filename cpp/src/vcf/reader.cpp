// All content Copyright (C) 2018 Genomics plc
#include <boost/algorithm/string/predicate.hpp>
#include <boost/regex.hpp>
#include <utils/exceptions.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <variant/type/variant.hpp>
#include <boost/lexical_cast.hpp>

#include "reader.hpp"

namespace wecall
{
namespace vcf
{
    Info parseVCFInfo( const std::string & raw_info )
    {
        Info parsed_info;

        std::vector< std::string > fields;
        boost::split( fields, raw_info, boost::is_any_of( ";" ) );

        for ( std::string & field : fields )
        {
            std::vector< std::string > parts;
            boost::split( parts, field, boost::is_any_of( "=" ) );

            if ( parts.size() == 1 )
            {
                parsed_info.push_back( std::make_pair( parts.front(), std::vector< std::string >{} ) );
            }
            else if ( parts.size() == 2 )
            {
                std::vector< std::string > value;
                boost::split( value, parts.back(), boost::is_any_of( "," ) );

                const auto notEmpty = []( const std::string & item )
                {
                    return not item.empty();
                };

                std::vector< std::string > cleaned_value;
                std::copy_if( value.begin(), value.end(), std::back_inserter( cleaned_value ), notEmpty );

                parsed_info.push_back( std::make_pair( parts.front(), cleaned_value ) );
            }
            else
            {
                WECALL_LOG( DEBUG, "Invalid info field " << field );
            }
        }

        return parsed_info;
    }
}
}
