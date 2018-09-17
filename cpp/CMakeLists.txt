cmake_minimum_required(VERSION 2.8)
project(wecall)

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(VENDOR "${CMAKE_CURRENT_SOURCE_DIR}/../vendor")
set(ENV{PROJECT_HOME} "${CMAKE_CURRENT_SOURCE_DIR}/..")

set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

set(COMMON_FLAGS "-std=c++11 -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -DBOOST_ALL_DYN_LINK --pedantic -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${COMMON_FLAGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${COMMON_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${COMMON_FLAGS}")

message("CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG}")
message("CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE}")
message("CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

execute_process(
    COMMAND bash -c "git rev-parse HEAD || echo none"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE GIT_SHA1
)

execute_process(
        COMMAND date -Is
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_STRIP_TRAILING_WHITESPACE
        OUTPUT_VARIABLE BUILD_DATE
)

execute_process(
        COMMAND cat ../VERSION
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_STRIP_TRAILING_WHITESPACE
        OUTPUT_VARIABLE PRODUCT_VERSION
)

configure_file("${CMAKE_SOURCE_DIR}/src/version/version.cpp.template" "${CMAKE_BINARY_DIR}/src/version/version.cpp" @ONLY)

set(WECALL_SOURCES
    src/weCallBase.hpp
    src/weCallMapAndReduce.cpp
    src/weCallMapAndReduce.hpp
    src/weCallReduce.cpp
    src/weCallReduce.hpp
    src/alignment/align.cpp
    src/alignment/align.hpp
    src/alignment/aligner.cpp
    src/alignment/aligner.hpp
    src/alignment/alignScorer.hpp
    src/alignment/galign.cpp
    src/alignment/galign.hpp
    src/alignment/mmHelpers.hpp
        src/assembly/node.cpp
    src/assembly/node.hpp
    src/assembly/sequenceGraph.cpp
    src/assembly/sequenceGraph.hpp
    src/caller/diploid/diploid.cpp
    src/caller/diploid/diploid.hpp
    src/caller/diploid/diploidAnnotate.hpp
    src/caller/diploid/diploidAnnotate.cpp
    src/caller/diploid/genotypeUtils.hpp
    src/caller/diploid/genotypeUtils.cpp
    src/caller/diploid/readSupportAccountant.cpp
    src/caller/diploid/readSupportAccountant.hpp
    src/caller/diploid/referenceCalling.cpp
    src/caller/diploid/referenceCalling.hpp
    src/caller/diploid/variantQualityCalculator.cpp
    src/caller/diploid/variantQualityCalculator.hpp
    src/caller/alignPhasing.cpp
    src/caller/alignPhasing.hpp
    src/caller/annotation.cpp
    src/caller/annotation.hpp
    src/caller/callSet.cpp
    src/caller/callSet.hpp
    src/caller/candidateVariantBank.cpp
    src/caller/candidateVariantBank.hpp
    src/caller/job.cpp
    src/caller/job.hpp
    src/caller/jobReduce.cpp
    src/caller/jobReduce.hpp
    src/caller/haplotypeLikelihoods.cpp
    src/caller/haplotypeLikelihoods.hpp
    src/caller/metadata.hpp
    src/caller/params.cpp
    src/caller/params.hpp
    src/caller/region.cpp
    src/caller/region.hpp
    src/caller/regionUtils.cpp
    src/caller/regionUtils.hpp
    src/caller/typedAnnotation.hpp
    src/io/readfilters/baseQualityFilter.cpp
    src/io/readfilters/baseQualityFilter.hpp
    src/io/readfilters/booleanFilter.cpp
    src/io/readfilters/booleanFilter.hpp
    src/io/readfilters/mapQualityFilter.cpp
    src/io/readfilters/mapQualityFilter.hpp
    src/io/readfilters/rangeFilter.hpp
    src/io/readfilters/readFilter.cpp
    src/io/readfilters/readFilter.hpp
    src/io/readfilters/readFilterAndTrimmer.cpp
    src/io/readfilters/readFilterAndTrimmer.hpp
    src/io/readfilters/shortReadFilter.cpp
    src/io/readfilters/shortReadFilter.hpp
    src/io/bamFile.cpp
    src/io/bamFile.hpp
    src/io/bamFileIterator.cpp
    src/io/bamFileIterator.hpp
    src/io/bedFile.cpp
    src/io/bedFile.hpp
    src/io/fastaFile.cpp
    src/io/fastaFile.hpp
    src/io/pysam.cpp
    src/io/pysam.hpp
    src/io/read.cpp
    src/io/read.hpp
    src/io/readIntervalTree.hpp
    src/io/readRange.cpp
    src/io/readRange.hpp
    src/io/readDataReader.cpp
    src/io/readDataReader.hpp
    src/io/readDataSet.cpp
    src/io/readDataSet.hpp
    src/io/readSummaries.cpp
    src/io/readSummaries.hpp
    src/io/readUtils.hpp
    src/io/readUtils.cpp
    src/io/tabixFile.hpp
    src/io/tabixFile.cpp
    src/io/tabixVCFFile.cpp
    src/io/tabixVCFFile.hpp
    src/io/vcfWriter.cpp
    src/io/vcfWriter.hpp
    src/mapping/hashMapper.cpp
    src/mapping/hashMapper.hpp
    src/readrecalibration/commonTypes.hpp
    src/readrecalibration/errorCorrectionParameters.hpp
    src/readrecalibration/intermediateOutputWriter.hpp
    src/readrecalibration/readRecalibration.cpp
    src/readrecalibration/readRecalibration.hpp
    src/readrecalibration/intermediateOutputWriter.cpp
    src/readrecalibration/kmerDistribution.hpp
    src/readrecalibration/kmerDistribution.cpp
    src/readrecalibration/readDataForErrorPosterior.cpp
    src/readrecalibration/readDataForErrorPosterior.hpp
    src/readrecalibration/siteKmerDistribution.cpp
    src/readrecalibration/siteKmerDistribution.hpp
    src/readrecalibration/siteReadDataForErrorPosterior.cpp
    src/readrecalibration/siteReadDataForErrorPosterior.hpp
    src/stats/functions.cpp
    src/stats/functions.hpp
    src/stats/models.cpp
    src/stats/models.hpp
    src/utils/referenceSequence.cpp
    src/utils/referenceSequence.hpp
    src/utils/bestScoreSelector.cpp
    src/utils/bestScoreSelector.hpp
    src/utils/combinationGenerator.hpp
    src/utils/combinations.hpp
    src/utils/date.hpp
    src/utils/indexedProduct.hpp
    src/utils/matrix.cpp
    src/utils/matrix.hpp
    src/utils/identity.hpp
    src/utils/exceptions.hpp
    src/utils/flatten.hpp
    src/utils/partition.hpp
    src/utils/intervalTree.hpp
    src/utils/interval.cpp
    src/utils/interval.hpp
    src/utils/logging.hpp
    src/utils/logging.cpp
    src/utils/median.hpp
    src/utils/multinomialCoefficients.cpp
    src/utils/multinomialCoefficients.hpp
    src/utils/NeedlemanWunsch.cpp
    src/utils/NeedlemanWunsch.hpp
    src/utils/sequence.cpp
    src/utils/sequence.hpp
    src/utils/timer.hpp
    src/utils/timer.cpp
    src/utils/write.hpp
    src/varfilters/filter.cpp
    src/varfilters/filter.hpp
    src/varfilters/variantSoftFilterBank.cpp
    src/varfilters/variantSoftFilterBank.hpp
    src/alignment/cigar.cpp
    src/alignment/cigar.hpp
    src/alignment/cigarItems.cpp
    src/alignment/cigarItems.hpp
    src/variant/breakpointVariantGenerator.cpp
    src/variant/breakpointVariantGenerator.hpp
    src/variant/clustering.cpp
    src/variant/clustering.hpp
    src/variant/genotype.cpp
    src/variant/genotype.hpp
    src/variant/haplotype.cpp
    src/variant/haplotype.hpp
    src/variant/haplotypeGenerator.hpp
    src/variant/haplotypeGenerator.cpp
    src/variant/haplotypeRanker.hpp
    src/variant/haplotypeRanker.cpp
    src/variant/type/breakpoint.cpp
    src/variant/type/breakpoint.hpp
    src/variant/type/variant.cpp
    src/variant/type/variant.hpp
    src/variant/variantCombinations.cpp
    src/variant/variantCombinations.hpp
    src/variant/variantContainer.cpp
    src/variant/variantContainer.hpp
    src/variant/variantFilter.cpp
    src/variant/variantFilter.hpp
    src/variant/variantGenerator.cpp
    src/variant/variantGenerator.hpp
    src/variant/variantGenerationData.cpp
    src/variant/variantGenerationData.hpp
    src/variant/variantNormalizer.cpp
    src/variant/variantNormalizer.hpp
    src/variant/snpFinder.hpp
    src/variant/snpFinder.cpp
    src/vcf/field.cpp
    src/vcf/field.hpp
    src/vcf/filterDescription.cpp
    src/vcf/filterDescription.hpp
    src/vcf/header.cpp
    src/vcf/header.hpp
    src/vcf/reader.cpp
    src/vcf/reader.hpp
    src/vcf/record.cpp
    src/vcf/record.hpp
    src/common.hpp
    src/version/version.hpp
    src/version/version.cpp
)

set(UNITTEST_SOURCES
    test/unittest/alignment/testGAlign.cpp
        test/unittest/assembly/testNode.cpp
    test/unittest/assembly/testSequenceGraph.cpp
        test/unittest/caller/testAlignPhasing.cpp
    test/unittest/caller/testCandidateVariantBank.cpp
    test/unittest/caller/testParams.cpp
    test/unittest/caller/testRegion.cpp
    test/unittest/caller/testRegionUtils.cpp
    test/unittest/mapping/testHashMapper.cpp
    test/unittest/readrecalibration/testCommonTypes.cpp
    test/unittest/readrecalibration/testReadRecalibration.cpp
    test/unittest/readrecalibration/testSiteKmerDistribution.cpp
    test/unittest/readFilters/testReadFilterAndTrimmer.cpp
    test/unittest/stats/testFunctions.cpp
    test/unittest/alignment/testCigar.cpp
    test/unittest/alignment/testCigarItems.cpp
    test/unittest/variant/testBreakpointVariantGenerator.cpp
    test/unittest/variant/testBreakpoint.cpp
    test/unittest/variant/testHaplotype.cpp
    test/unittest/variant/testHaplotypeVector.cpp
    test/unittest/variant/testClustering.cpp
    test/unittest/variant/testDeletion.cpp
    test/unittest/variant/testGenotype.cpp
    test/unittest/variant/testGenotypeVector.cpp
    test/unittest/variant/testGenotypeVectorForPloidy1.cpp
    test/unittest/variant/testGenotypeVectorForPloidy2.cpp
    test/unittest/variant/testGenotypeVectorForPloidy3.cpp
    test/unittest/variant/testInsertion.cpp
    test/unittest/variant/testVariant.cpp
    test/unittest/variant/testMnp.cpp
    test/unittest/variant/testSnp.cpp
    test/unittest/variant/testHaplotypeGeneration.cpp
    test/unittest/variant/testHaplotypeRanker.cpp
    test/unittest/variant/testVariantCombinations.cpp
    test/unittest/variant/testVariantContainer.cpp
    test/unittest/variant/testVariantGenerator.cpp
    test/unittest/variant/testSnpFinder.cpp
    test/unittest/variant/testVariantNormalizer.cpp
    test/unittest/varfilters/testVariantSoftFilterBank.cpp
    test/unittest/utils/testBestScoreSelector.cpp
    test/unittest/utils/testReferenceSequence.cpp
    test/unittest/utils/testFactorial.cpp
    test/unittest/utils/testFlatten.cpp
    test/unittest/utils/testInterval.cpp
    test/unittest/utils/testMatrix.cpp
    test/unittest/utils/testMedian.cpp
    test/unittest/utils/testMultinomialCoefficients.cpp
    test/unittest/utils/testNeedlemanWunsch.cpp
    test/unittest/utils/testNWPenalties.cpp
    test/unittest/utils/testNWVariant.cpp
    test/unittest/utils/testPartition.cpp
    test/unittest/utils/testSequence.cpp
    test/unittest/utils/testWrite.cpp
    test/unittest/vcf/testField.cpp
    test/unittest/vcf/testRecord.cpp
    test/unittest/vcf/testReader.cpp
    test/unittest/vcf/VCFTestUtils.cpp
    test/unittest/vcf/VCFTestUtils.hpp
    test/unittest/unittest.cpp
    test/unittest/variant/testVariantFilter.cpp
#    test/unittest/alignment/testMMHelpers.cpp
    )

set(IOTEST_SOURCES
        test/ioTest/caller/testRegionUtils.cpp
        test/ioTest/io/ioFixture.hpp
        test/ioTest/io/testBedFile.cpp
        test/ioTest/io/testBuildRefCall.cpp
        test/ioTest/io/testReadDataset.cpp
        test/ioTest/io/testFastaFile.cpp
        test/ioTest/io/testRead.cpp
        test/ioTest/io/testReadRange.cpp
        test/ioTest/io/testReadIntervalTree.cpp
        test/ioTest/io/testReadUtils.cpp
        test/ioTest/io/testReadSummaries.cpp
        test/ioTest/io/testTabixFile.cpp
        test/ioTest/io/testVCFWriter.cpp
        test/ioTest/utils/environment.cpp
        test/ioTest/utils/environment.hpp
        test/ioTest/readrecalibration/testKmerDistrubution.cpp
        test/ioTest/ioTest.cpp)

include_directories(src)
include_directories(test)

include_directories(SYSTEM "${VENDOR}")

##################################
# Find boost
set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost COMPONENTS regex system filesystem program_options iostreams timer chrono date_time log log_setup thread unit_test_framework REQUIRED)
if ("${CMAKE_BUILD_TYPE}" STREQUAL "DEBUG")
    set(BOOST_SHARED_REGEX_LIBRARY ${Boost_REGEX_LIBRARY_DEBUG})
    set(BOOST_SHARED_SYSTEM_LIBRARY ${Boost_SYSTEM_LIBRARY_DEBUG})
    set(BOOST_SHARED_FILESYSTEM_LIBRARY ${Boost_FILESYSTEM_LIBRARY_DEBUG})
    set(BOOST_SHARED_PROGRAM_OPTIONS_LIBRARY ${Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG})
    set(BOOST_SHARED_IOSTREAMS_LIBRARY ${Boost_IOSTREAMS_LIBRARY_DEBUG})
    set(BOOST_SHARED_TIMER_LIBRARY ${Boost_TIMER_LIBRARY_DEBUG})
    set(BOOST_SHARED_CHRONO_LIBRARY ${Boost_CHRONO_LIBRARY_DEBUG})
    set(BOOST_SHARED_DATE_TIME_LIBRARY ${Boost_DATE_TIME_LIBRARY_DEBUG})
    set(BOOST_SHARED_LOG_LIBRARY ${Boost_LOG_LIBRARY_DEBUG})
    set(BOOST_SHARED_LOG_SETUP_LIBRARY ${Boost_LOG_SETUP_LIBRARY_DEBUG})
    set(BOOST_SHARED_THREAD_LIBRARY ${Boost_THREAD_LIBRARY_DEBUG})
    set(BOOST_SHARED_UNIT_TEST_FRAMEWORK_LIBRARY ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG})
else()
    set(BOOST_SHARED_REGEX_LIBRARY ${Boost_REGEX_LIBRARY_RELEASE})
    set(BOOST_SHARED_SYSTEM_LIBRARY ${Boost_SYSTEM_LIBRARY_RELEASE})
    set(BOOST_SHARED_FILESYSTEM_LIBRARY ${Boost_FILESYSTEM_LIBRARY_RELEASE})
    set(BOOST_SHARED_PROGRAM_OPTIONS_LIBRARY ${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE})
    set(BOOST_SHARED_IOSTREAMS_LIBRARY ${Boost_IOSTREAMS_LIBRARY_RELEASE})
    set(BOOST_SHARED_TIMER_LIBRARY ${Boost_TIMER_LIBRARY_RELEASE})
    set(BOOST_SHARED_CHRONO_LIBRARY ${Boost_CHRONO_LIBRARY_RELEASE})
    set(BOOST_SHARED_DATE_TIME_LIBRARY ${Boost_DATE_TIME_LIBRARY_RELEASE})
    set(BOOST_SHARED_LOG_LIBRARY ${Boost_LOG_LIBRARY_RELEASE})
    set(BOOST_SHARED_LOG_SETUP_LIBRARY ${Boost_LOG_SETUP_LIBRARY_RELEASE})
    set(BOOST_SHARED_THREAD_LIBRARY ${Boost_THREAD_LIBRARY_RELEASE})
    set(BOOST_SHARED_UNIT_TEST_FRAMEWORK_LIBRARY ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE})
endif()

include(FindZLIB)

SET(BOOST_LIBRARIES
        ${BOOST_SHARED_REGEX_LIBRARY}
        ${BOOST_SHARED_SYSTEM_LIBRARY}
        ${BOOST_SHARED_FILESYSTEM_LIBRARY}
        ${BOOST_SHARED_PROGRAM_OPTIONS_LIBRARY}
        ${BOOST_SHARED_IOSTREAMS_LIBRARY}
        ${BOOST_SHARED_TIMER_LIBRARY}
        ${BOOST_SHARED_CHRONO_LIBRARY}
        ${BOOST_SHARED_DATE_TIME_LIBRARY}
        ${BOOST_SHARED_LOG_LIBRARY}
        ${BOOST_SHARED_LOG_SETUP_LIBRARY}
        ${BOOST_SHARED_THREAD_LIBRARY}
        ${BOOST_SHARED_UNIT_TEST_FRAMEWORK_LIBRARY}
        )
SET(BOOST_DEB_VERSION "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
message("BOOST_LIBRARIES ${BOOST_LIBRARIES}")

set(LIBS
    pthread
    bam
    tabix
    ${ZLIB_LIBRARIES}
)
#set(COMMON_LINK_DIRECTORIES "${PROJECT_ENV}/lib" "${DEPENDENCIES}/libboost/lib" "${DEPENDENCIES}/samtools/lib" "${DEPENDENCIES}/tabix/lib")
set(COMMON_LINK_DIRECTORIES 
  "${VENDOR}/tabix" 
  "${VENDOR}/samtools"
)

foreach(library ${LIBS})
    find_library(${library} "${COMMON_LINK_DIRECTORIES}")
endforeach(library)

set(LINK_DIRECTORIES
    "${COMMON_LINK_DIRECTORIES}"
)
link_directories(${LINK_DIRECTORIES})

add_executable(weCall ${WECALL_SOURCES} src/weCall.cpp)
target_link_libraries(weCall ${BOOST_LIBRARIES} ${LIBS})

add_executable(unittest ${WECALL_SOURCES} ${UNITTEST_SOURCES})
target_link_libraries(unittest ${BOOST_LIBRARIES} ${LIBS})

add_executable(iotest ${WECALL_SOURCES} ${IOTEST_SOURCES})
target_link_libraries(iotest ${BOOST_LIBRARIES} ${LIBS})

install(TARGETS weCall unittest iotest RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
