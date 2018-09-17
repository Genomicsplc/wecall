// All content Copyright (C) 2018 Genomics plc
#ifndef CALLER_CALLSET_HPP
#define CALLER_CALLSET_HPP

#include <limits>

#include "common.hpp"

#include "variant/type/variant.hpp"
#include "caller/annotation.hpp"

namespace echidna
{
namespace caller
{
    struct Call : public Annotated
    {
        enum Type
        {
            VAR,
            REF,
            UNKNOWN
        };

        Type type;
        variant::varPtr_t var;
        utils::Interval interval;
        double qual;
        std::set< std::string > filters;
        std::set< std::string > varIds;

        struct Sample : public Annotated
        {
            std::vector< Type > genotypeCalls;

            bool hasVar() const
            {
                return std::find( genotypeCalls.begin(), genotypeCalls.end(), VAR ) != genotypeCalls.end();
            }
        };
        std::vector< Sample > samples;

        Call( const variant::varPtr_t & theVar,
              const utils::Interval theInterval,
              double theVarQ,
              std::size_t nSamples,
              const std::vector< std::vector< Type > > & genoCalls )
            : type( theVar != nullptr ? VAR : REF ),
              var( theVar ),
              interval( theInterval ),
              qual( theVarQ ),
              samples( nSamples )
        {
            for ( std::size_t i = 0; i < nSamples; ++i )
            {
                samples[i].genotypeCalls = genoCalls[i];
            }
        }

        //        Call( const Call & other ) = default;
        //        Call( Call && other ) = default;

        bool isRefCall() const { return type == REF; }
    };

    struct CallComp
    {
        bool operator()( const Call & x, const Call & y ) const { return variant::varPtrComp()( x.var, y.var ); }
    };

    using callVector_t = std::vector< Call >;
}

// Injected into project namespace - to be used consistently across the project

using genoCall_t = std::vector< caller::Call::Type >;
using genoCalls_t = std::vector< std::vector< caller::Call::Type > >;

namespace caller
{
    bool variantCalled( const genoCalls_t & genoCalls );
    int64_t zeroIndexedVCFPosition( const callVector_t & genoCall );
}

using callIt_t = caller::callVector_t::const_iterator;
}

#endif
