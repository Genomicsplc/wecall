// All content Copyright (C) 2018 Genomics plc
#ifndef TYPED_ANNOTATION_HPP
#define TYPED_ANNOTATION_HPP

#include <string>
#include <memory>
#include <vector>
#include <cassert>
#include <iostream>
#include <cmath>
#include <type_traits>

#include "caller/annotation.hpp"

namespace echidna
{
namespace caller
{
    /// Dummy function defined so I can use std::to_string with string types
    std::string to_string( const std::string & theString ) { return theString; }

    ///---------------------------------------------------------------------------------------------

    template < typename T >
    class TypedAnnotation : public Annotation
    {
    public:
        TypedAnnotation( const std::string & id, const T value ) : m_id( id ), m_value( value ) {}

        virtual const std::string & getID() const { return m_id; }

        virtual std::string getValue() const { return to_string( m_value ); }

    private:
        std::string m_id;
        T m_value;
    };

    //-------------------------------------------------------------------------------------------------
}
}

#endif
