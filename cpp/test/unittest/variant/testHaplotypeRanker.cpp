#include "mapping/hashMapper.hpp"
#include "utils/exceptions.hpp"
#include "variant/type/variant.hpp"
#include "utils/referenceSequence.hpp"
#include "caller/region.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <string>
#include "io/read.hpp"
#include "variant/haplotypeRanker.hpp"
#include <boost/algorithm/string.hpp>
