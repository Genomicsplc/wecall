// All content Copyright (C) 2018 Genomics plc
#ifndef BASEQUAL_FILTER_HPP
#define BASEQUAL_FILTER_HPP

#include "io/readfilters/readFilter.hpp"

namespace echidna
{
namespace io
{
    /// Filter which checks the base qualities of the specified read.
    class BaseQualFilter : public ReadFilter
    {
    public:
        /// Constructor.
        ///
        /// @param threshold The minimum acceptable base quality value
        /// @param minBases The minimum number of bases that must have Qual > threshold
        BaseQualFilter( int threshold, int minBases ) : m_qualThreshold( threshold ), m_minBases( minBases ) {}

        /// Destructor
        virtual ~BaseQualFilter(){};

        /// Convert to string representation
        std::string toString() const override;

    private:
        virtual bool passesFilter_impl( const io::Read & theRead ) override;
        int m_qualThreshold;  /// Min value for base quality.
        int m_minBases;       /// This many bases must have higher qual than threshold
    };
}
}

#endif
