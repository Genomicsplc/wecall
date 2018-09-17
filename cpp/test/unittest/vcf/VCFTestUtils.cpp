#include "VCFTestUtils.hpp"
#include <boost/test/unit_test.hpp>
#include <typeinfo>
#include <typeindex>

void check_equal_info( const echidna::vcf::Info & expected, const echidna::vcf::Info & parsed )
{
    std::ostringstream info_size_message;
    info_size_message << "Different length info fields (expected " << expected.size() << " elements: [";
    for ( auto & item : expected )
    {
        info_size_message << item.first << ", ";
    }
    info_size_message << "]; got " << parsed.size() << " elements: [";
    for ( auto & item : parsed )
    {
        info_size_message << item.first << ", ";
    }
    info_size_message << "]).";
    BOOST_CHECK_MESSAGE( expected.size() == parsed.size(), info_size_message.str().c_str() );
    for ( unsigned int i = 0; i != std::min( expected.size(), parsed.size() ); ++i )
    {
        std::ostringstream item_key_message;
        item_key_message << "Different keys (expected '" << expected[i].first << "'; got '" << parsed[i].first
                         << "') for i = " << i;
        BOOST_CHECK_MESSAGE( expected[i].first == parsed[i].first, item_key_message.str().c_str() );
        std::ostringstream inner_length_message;
        inner_length_message << "Inner vectors of different size (expected '" << expected[i].second.size() << "'; got '"
                             << parsed[i].second.size() << "') for i = " << i;
        BOOST_CHECK_MESSAGE( expected[i].second.size() == parsed[i].second.size(), inner_length_message.str().c_str() );
        for ( unsigned int j = 0; j != std::min( expected[i].second.size(), parsed[i].second.size() ); ++j )
        {
            std::ostringstream inner_element_message;
            inner_element_message << "Elements don't match (expected '" << expected[i].second[j] << "'; got '"
                                  << parsed[i].second[j] << "') for i, j = " << i << ", " << j;
            BOOST_CHECK_MESSAGE( expected[i].second[j] == parsed[i].second[j], inner_element_message.str().c_str() );
        }
    }
}

//-------------------------------------------------------------------------------------------------

boost::test_tools::predicate_result checkVariantInVector( const std::vector< varPtr_t > & variants,
                                                          const varPtr_t & testVariant )
{
    auto comparison_func = [&testVariant]( const varPtr_t & v )
    {

        return std::type_index( typeid( *testVariant ) ) == std::type_index( typeid( *v ) ) and *v == *testVariant;
    };

    const auto it = std::find_if( variants.cbegin(), variants.cend(), comparison_func );
    boost::test_tools::predicate_result res( false );
    res.message() << "\nVariant " << testVariant->toString().c_str() << " not found in the set:\n{\n";
    for ( auto & var : variants )
    {
        res.message() << "\t" << var->toString().c_str() << "\n";
    }
    res.message() << "}";

    if ( it == variants.end() )
    {
        return res;
    }
    // else
    return true;
}
