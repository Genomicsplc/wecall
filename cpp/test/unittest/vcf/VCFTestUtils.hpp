// All content Copyright (C) 2018 Genomics plc
#ifndef VCF_VCFTESTUTILS_HPP
#define VCF_VCFTESTUTILS_HPP

#include <boost/test/predicate_result.hpp>
#include "vcf/record.hpp"

using wecall::variant::varPtr_t;

boost::test_tools::predicate_result checkVariantInVector( const std::vector< varPtr_t > & variants,
                                                          const varPtr_t & testVariant );

#endif
