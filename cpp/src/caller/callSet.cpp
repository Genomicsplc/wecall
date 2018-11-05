// All content Copyright (C) 2018 Genomics plc
#include "caller/callSet.hpp"

namespace wecall
{
namespace caller
{
    bool variantCalled( const genoCalls_t & genoCalls )
    {
        bool called = false;
        for ( const auto & genoCall : genoCalls )
        {
            for ( const auto & call : genoCall )
            {
                if ( call == Call::VAR )
                {
                    called = true;
                }
            }
        }
        return called;
    }

    int64_t zeroIndexedVCFPosition( const callVector_t & calls )
    {
        auto position = std::numeric_limits< int64_t >::max();
        auto variantExists = false;
        for ( const auto & call : calls )
        {
            if ( not call.isRefCall() )
            {
                variantExists = true;
                position = std::min( call.var->zeroIndexedVcfPosition(), position );
            }
        }
        if ( variantExists )
        {
            return position;
        }
        else
        {
            return -1;
        }
    }

    //-----------------------------------------------------------------------------------------
}
}
