// All content Copyright (C) 2018 Genomics plc
#ifndef ANNOTATION_HPP
#define ANNOTATION_HPP

#include "common.hpp"
#include "utils/write.hpp"
#include "utils/exceptions.hpp"

#include <iostream>
#include <memory>
#include <vector>
#include <vcf/field.hpp>

namespace echidna
{
namespace caller
{
    ///---------------------------------------------------------------------------------------------

    struct Annotation
    {
        enum Type
        {
            COUNT,
            PHRED_VAL,
            LOG_VAL,
            DOUBLE_VAL,
            ID,
            FLAG
        };

        template < typename T >
        struct Def
        {
            const std::string id;
            const Type type;

            Def( const std::string & inId, Type inType ) : id( inId ), type( inType ) {}
        };

        virtual ~Annotation() {}
        virtual std::string getID() const = 0;
        virtual std::vector< std::string > getValues() const = 0;

        // TODO - These should be tied in with the list of IDs in vcfwriter.cpp

        // Call annotations...
        static const Def< int64_t > BEG, END, LEN, DP, DPR, DPF, VC, VCR, VCF;
        static const Def< phred_t > PP, BR;
        static const Def< double > ABPV, SBPV, MQ, QD;

        // Genotype call annotations...
        static const Def< int64_t > PS, MIN_DP, FORMAT_DP;
        static const Def< phred_t > GQ, PQ;
        static const Def< std::vector< phred_t > > PL;
        static const Def< std::vector< int64_t > > AD;
        static const Def< std::vector< double > > VAF;

        template < typename T >
        static Def< T > define( std::string id, Type type )
        {
            return Def< T >( id, type );
        }
    };

    template < typename T >
    struct TypedAnnotation : Annotation
    {
        TypedAnnotation( const Def< T > & inDef, T inData ) : def( inDef ), data( inData ) {}

        std::string getID() const { return def.id; }
        std::vector< std::string > getValues() const;

        const Def< T > def;
        T data;
    };

    std::string serialise_double( const double data, const Annotation::Type type );
    std::string serialise_phred( const phred_t data );
    std::string serialise_integral_type( const int64_t data );
}

// Injected into project namespace - to be used consistently across the project.
using AnnotationPtr_t = std::shared_ptr< caller::Annotation >;
using Annotation = caller::Annotation;

namespace caller
{
    //-------------------------------------------------------------------------------------------------

    class Annotated
    {
    private:
        std::vector< AnnotationPtr_t > m_annotations;

    public:
        template < typename T >
        void addAnnotation( const Annotation::Def< T > & def, T data )
        {
            m_annotations.emplace_back( std::make_shared< TypedAnnotation< T > >( def, data ) );
        }

        template < typename T >
        T & getAnnotation( const Annotation::Def< T > & def ) const
        {
            for ( const AnnotationPtr_t & annotation : m_annotations )
            {
                if ( annotation->getID() == def.id )
                {
                    return ( std::dynamic_pointer_cast< TypedAnnotation< T > >( annotation ) )->data;
                }
            }

            // Oops! Cannot find annotation!
            throw utils::echidna_exception( "Cannot find annotation ID: " + def.id );
        }

        const std::vector< AnnotationPtr_t > & getAnnotations() const { return m_annotations; }
    };
}

//-------------------------------------------------------------------------------------------------
}

#endif
