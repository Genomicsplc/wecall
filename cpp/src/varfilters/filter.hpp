// All content Copyright (C) 2018 Genomics plc
#ifndef VARFILTERS_FILTER_HPP
#define VARFILTERS_FILTER_HPP

#include "common.hpp"
#include "vcf/field.hpp"
#include "vcf/filterDescription.hpp"
#include "caller/callSet.hpp"
#include "variant/type/variant.hpp"

#include <string>
#include <iostream>
#include <unordered_set>

namespace wecall
{
namespace varfilters
{
    /// The base class for soft variant output filters.
    class Filter
    {
    public:
        /// Basic constructor from constituent data.
        ///
        /// @param id Filter ID.
        /// @param description Filter description.
        Filter( const std::string & id, const std::string & description );

        /// Destructor declared virtual as may be sub-classed.
        virtual ~Filter() {}

        /// Returns the filter ID
        ///
        /// @return Filter ID.
        std::string getID() { return m_id; }

        /// Checks if a call is caught by (fails) this filter
        ///
        /// @param call Call to be tested
        /// @return True if the call is caught by (fails) this filter
        virtual bool catches( const caller::Call & call );

        /// Returns the FilterDesc object for this filter
        vcf::FilterDesc getVCFFilterDesc();

        /// Writes the filter to the output stream in.
        ///
        /// @param out Output stream
        /// @param filter Filter to be output
        /// @return Output stream
        friend std::ostream & operator<<( std::ostream & out, const Filter & filter );

    private:
        std::string m_id;
        std::string m_description;
    };

    // Filter pointers...
    using FilterPtr_t = std::shared_ptr< Filter >;

    // Comparison function used to sort shared pointers to filter objects
    struct FilterPtrComp
    {
        bool operator()( const FilterPtr_t & x, const FilterPtr_t & y ) const;
    };

    /// An output filter that tests allele bias.
    class ABFilter : public Filter
    {
    public:
        /// Basic constructor
        ///
        /// @param thresholdP probability threshold below which allele bias will be caught
        explicit ABFilter( double thresholdP );

        virtual ~ABFilter() {}
        virtual bool catches( const caller::Call & call ) override;

    private:
        double m_thresholdP;
    };

    /// An output filter that tests strand bias.
    class SBFilter : public Filter
    {
    public:
        /// Basic constructor
        ///
        /// @param thresholdP probability threshold below which strand bias will be caught
        explicit SBFilter( double thresholdP );

        virtual ~SBFilter() {}
        virtual bool catches( const caller::Call & call ) override;

    private:
        double m_thresholdP;
    };

    /// An output filter that tests allele + strand bias.
    class ABPlusSBFilter : public Filter
    {
    public:
        explicit ABPlusSBFilter( double thresholdP );

        virtual ~ABPlusSBFilter() {}
        virtual bool catches( const caller::Call & call ) override;

    private:
        double m_thresholdP;
    };

    class MQFilter : public Filter
    {
    public:
        /// Basic constructor
        ///
        explicit MQFilter( double thresholdP );

        virtual ~MQFilter() {}
        virtual bool catches( const caller::Call & call ) override;

    private:
        double m_thresholdP;
    };

    class QDFilter : public Filter
    {
    public:
        /// Basic constructor
        ///
        explicit QDFilter( double SNPThresholdP, double INDELThresholdP );

        virtual ~QDFilter() {}
        virtual bool catches( const caller::Call & call ) override;

    private:
        double m_SNPThresholdP;
        double m_INDELThresholdP;
    };

    class BRFilter : public Filter
    {
    public:
        /// Basic constructor
        ///
        explicit BRFilter( phred_t thresholdPhred );

        virtual ~BRFilter() {}
        virtual bool catches( const caller::Call & call ) override;

    private:
        phred_t m_thresholdPhred;
    };

    class LQFilter : public Filter
    {
    public:
        /// Basic constructor
        ///
        explicit LQFilter( phred_t thresholdPhred );

        virtual ~LQFilter() {}
        virtual bool catches( const caller::Call & call ) override;

    private:
        phred_t m_thresholdPhred;
    };
}
}

#endif
