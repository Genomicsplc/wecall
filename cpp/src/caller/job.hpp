// All content Copyright (C) 2018 Genomics plc
#ifndef CALLER_JOB_HPP
#define CALLER_JOB_HPP

#include "caller/diploid/diploid.hpp"
#include "common.hpp"
#include "caller/region.hpp"
#include "caller/params.hpp"
#include "caller/candidateVariantBank.hpp"
#include "io/readDataReader.hpp"
#include "io/readRange.hpp"
#include "io/fastaFile.hpp"
#include "io/vcfWriter.hpp"
#include "io/tabixVCFFile.hpp"
#include "readrecalibration/intermediateOutputWriter.hpp"
#include "variant/clustering.hpp"
#include "variant/type/variant.hpp"
#include "varfilters/variantSoftFilterBank.hpp"

namespace echidna
{
namespace caller
{
    class Job
    {
    public:
        Job( caller::params::Application applicationParams,
             caller::params::Data dataParams,
             caller::params::System systemParams,
             caller::params::PrivateSystem privateSystemParams,
             caller::params::Filters filterParams,
             caller::params::PrivateCalling privateCallingParams,
             caller::params::Calling callingParams,
             caller::params::PrivateData privateDataParams );

        void process();

    private:
        /// Second tier of processing - job is split into manageable blocks
        /// within process() and each is processed in turn. Within a block, reads
        /// are filtered and candidate variants generated from them. These are then
        /// split into clusters from which variants are called.
        int64_t processBlock( io::readDataset_t readDataset,
                              bool lastBlock,
                              const std::vector< std::size_t > & ploidyPerSample );

        callVector_t processBigCluster( const variant::VariantCluster & cluster,
                                        io::readDataset_t readDataset,
                                        const utils::referenceSequencePtr_t & referenceSequence,
                                        const std::vector< std::size_t > & ploidyPerSample );

        callVector_t processCluster( const variant::VariantCluster & cluster,
                                     const io::perSampleRegionsReads_t & regionReads,
                                     const io::perSampleRegionsReads_t & allReads,
                                     const utils::referenceSequencePtr_t & blockReferenceSequence,
                                     const std::vector< std::size_t > & ploidyPerSample );

        utils::referenceSequencePtr_t getReferenceForCluster(
            const variant::VariantCluster & cluster,
            const io::perSampleRegionsReads_t & regionReads,
            const utils::referenceSequencePtr_t & blockReferenceSequence ) const;

        void recalibrateBaseQualities( const io::perSampleRegionsReads_t & readRangesPerSample,
                                       const caller::Region & region ) const;

        std::vector< variant::VariantCluster > generateVariantClustersInBlock(
            const caller::Region & blockRegion,
            io::perSampleRegionsReads_t allReads,
            utils::referenceSequencePtr_t referenceSequence );

        void writeCallsForCluster( const callVector_t calls,
                                   const variant::VariantCluster cluster,
                                   const int64_t refStart,
                                   io::readDataset_t readDataset,
                                   const std::vector< std::size_t > & ploidyPerSample );

        void callReference( const std::string & contig,
                            int64_t start,
                            int64_t end,
                            io::readDataset_t readDataset,
                            const std::vector< std::size_t > & ploidyPerSample );

        callVector_t filterOutputCalls( const std::string & contig, const callVector_t & calls ) const;

        const caller::params::Application m_applicationParams;
        const caller::params::Data m_dataParams;
        const caller::params::System m_systemParams;
        const caller::params::PrivateSystem m_privateSystemParams;
        const caller::params::Filters m_filterParams;
        const caller::params::PrivateCalling m_privateCallingParams;
        const caller::params::Calling m_callingParams;
        const caller::params::PrivateData m_privateDataParams;

        std::shared_ptr< io::TabixVCFFile > m_candidateVCFFile;
        std::shared_ptr< io::TabixVCFFile > m_genotypeAllelesVCFFile;

        corrector::IntermediateOutputWriter m_intermediateOutputWriter;

        varfilters::VariantSoftFilterBank m_variantSoftFilterBank;
        io::ReadDataReader m_readDataReader;
        io::VCFWriter m_vcOut;
        io::FastaFile m_ref;

        const regions_t m_outputRegions;
        const regions_t m_callingRegions;

        const model::Model m_model;
        const model::VCFCallVectorBuilder m_variantCallBuilder;
    };
}
}

#endif
