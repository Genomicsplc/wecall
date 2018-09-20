#include "VCFTestUtils.hpp"
#include <boost/test/unit_test.hpp>
#include <typeinfo>
#include <typeindex>

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
