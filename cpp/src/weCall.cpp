// All content Copyright (C) 2018 Genomics plc
#include "weCallMapAndReduce.hpp"
#include "weCallReduce.hpp"

int main( const int argc, char * argv[] )
{
    if ( argc > 1 and strncmp( argv[1], "reduce", 7 ) == 0 )
    {
        // Shift command-line arguments up by one
        char ** newArgv = new char * [argc - 1];
        newArgv[0] = argv[0];

        for ( std::size_t i = 1; i < static_cast< std::size_t >( argc - 1 ); ++i )
        {
            newArgv[i] = argv[i + 1];
        }

        const auto rt = echidna::weCallReduce().processJob( argc - 1, newArgv );
        delete[] newArgv;
        return rt;
    }
    else
    {
        return echidna::weCallMapAndReduce().processJob( argc, argv );
    }
}
