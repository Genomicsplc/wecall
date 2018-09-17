# All content Copyright (C) 2018 Genomics plc
import logging
import re

from wecall.bamutils.bam_builder import BAMBuilder
from wecall.bedutils.bedwriter import BEDWriterContextManager
from wecall.wecall_utils.wecall_config_builder import WecallConfigBuilder
from wecall.wecall_utils.wecall_input_data_builder import WecallInputDataBuilder
from wecall_test_drivers.tool_runner import log_file, log_bam_file
import os
from wecall.bamutils.read_sequence import HIGH_QUALITY
from wecall.bamutils.sample_bank import SampleBank
from wecall.genomics.reference_chromosome import DEFAULT_CHROM
from wecall_test_drivers.ascii_wecall_runner import DEFAULT_SAMPLE_NAME
from wecall_test_drivers.variant_caller_wrapper import VariantCallerWrapper
from wecall_test_drivers.vcf_expectation import VCFExpectation


class SVCDriver(object):

    def __init__(self, test_case):
        self.__test_case = test_case
        self._sample_bank = {}
        self.__bam_data = {}
        self._config = {
            'noSimilarReadsFilter': False,
            'minCallQual': 2,
        }
        self._output_vcf_filename = None
        self.__log_filename = None
        self.__bam_filenames = None
        self.__reference_filename = None

    def with_ref_sequence(self, ref_sequence, pos_from=0, chrom=DEFAULT_CHROM):
        self._sample_bank[chrom] = SampleBank(ref_sequence, pos_from, chrom)
        return self

    def with_read(
            self,
            read,
            quality=None,
            n_fwd=None,
            n_rev=None,
            mapping_quality=HIGH_QUALITY,
            chrom=DEFAULT_CHROM,
            sample_name=DEFAULT_SAMPLE_NAME,
            read_id=None,
            read_flags=None,
            cigar_string=None,
            read_start=None,
            read_mate_start=None
    ):
        assert(chrom in self._sample_bank.keys())
        if sample_name not in self._sample_bank[chrom].sample_names:
            self._sample_bank[chrom].add_sample_name(sample_name)
        self._sample_bank[chrom][sample_name].add_sequence(
            read,
            n_fwd=n_fwd,
            n_rev=n_rev,
            mapping_quality=mapping_quality,
            quality_string=quality,
            read_id=read_id,
            read_flags=read_flags,
            cigar_string=cigar_string,
            read_start=read_start,
            read_mate_start=read_mate_start
        )
        return self

    def with_ploidy(self, ploidy):
        self._config['ploidy'] = ploidy
        return self

    def with_normalize_variant_calls(self, normalize_variant_calls_option):
        self._config['normalizeVariantCalls'] = normalize_variant_calls_option
        return self

    def with_simple_reads(self):
        # for adding reads without explicitly specified metadata
        raise NotImplementedError

    def with_ref_filename(self, filename):
        self.__reference_filename = filename
        return self

    def with_read_with_quality(
            self,
            read,
            quality,
            sample_name=None,
            n_fwd=None,
            n_rev=None):
        raise NotImplementedError()

    def with_bam_data(self, file_name, bam_data, with_read_group=True):
        self.__bam_data[file_name] = (bam_data, with_read_group)
        return self

    def with_output_vcf_filename(self, filename):
        self._output_vcf_filename = filename
        return self

    def with_min_call_qual(self, qual_phred):
        self._config['minCallQual'] = qual_phred
        return self

    def with_output_format(self, outputFormat):
        self._config['outputFormat'] = outputFormat
        return self

    def with_mem_limit(self, mem_limit):
        self._config['memLimit'] = mem_limit
        return self

    def with_all_variants(self, all_variants):
        self._config['allVariants'] = all_variants
        return self

    def with_candidate_variants_file(self, filename):
        self._config['candidateVariantsFile'] = filename
        return self

    def with_genotype_alleles(self, filename):
        self._config['genotypeAllelesFile'] = filename
        return self

    def with_allow_MNP_calls(self, condition):
        self._config['allowMNPCalls'] = condition
        return self

    def with_max_cluster_distance(self, distance):
        self._config['maxClusterDist'] = distance
        return self

    def with_min_cluster_distance(self, distance):
        self._config['minClusterDist'] = distance
        return self

    def with_output_phased_genotypes(self, condition):
        self._config["outputPhasedGenotypes"] = condition
        return self

    def with_allow_improper_pairs(self):
        self._config["allowImproperPairs"] = True
        return self

    def with_duplicates_filter(self, condition):
        self._config["duplicatesFilter"] = condition
        return self

    def with_no_similar_reads_filter(self, condition):
        self._config["noSimilarReadsFilter"] = condition
        return self

    def with_var_filters(self, *var_filters):
        self._config["varFilterIDs"] = ",".join(var_filters)
        return self

    def with_bad_reads_window_size(self, bad_reads_window_size):
        self._config['badReadsWindowSize'] = bad_reads_window_size
        return self

    def with_min_bad_reads_score(self, min_bad_reads_score):
        self._config['minBadReadsScore'] = min_bad_reads_score
        return self

    def with_min_snp_q_over_depth(self, min_quality_over_depth):
        self._config['minSNPQOverDepth'] = min_quality_over_depth
        return self

    def with_overwrite(self, status):
        self._config['overwrite'] = status
        return self

    def with_number_of_jobs(self, n_jobs):
        self._config['numberOfJobs'] = n_jobs
        return self

    def with_work_dir(self, location):
        self._config['workDir'] = location
        return self

    def with_min_indel_q_over_depth(self, min_quality_over_depth):
        self._config['minIndelQOverDepth'] = min_quality_over_depth
        return self

    def with_min_root_mean_square_mapping_q(self, min_rms_mapping_quality):
        self._config['minRMSMappingQ'] = min_rms_mapping_quality
        return self

    def with_strand_bias_p(self, p_value):
        self._config['minStrandBiasP'] = p_value
        return self

    def with_allele_plus_strand_bias_p(self, p_value):
        self._config['minAllelePlusStrandBiasP'] = p_value
        return self

    def with_read_mapping_filter_q(self, min_read_mapping_quality):
        self._config['readMappingFilterQ'] = min_read_mapping_quality
        return self

    def with_log_timings(self, condition):
        self._config['logTimings'] = condition
        return self

    def with_region_string(self, region_string):
        self._config['regions'] = region_string
        return self

    def with_region_padding(self, padding):
        self._config['regionPadding'] = padding
        return self

    def with_min_reads_per_var(self, value):
        self._config["minReadsPerVar"] = value
        return self

    def with_output_ref_calls(self, condition):
        self._config["outputRefCalls"] = condition
        return self

    def with_max_ref_call_size(self, size):
        self._config["maxRefCallSize"] = size
        return self

    def with_log_filename(self, filename):
        self.__log_filename = filename
        return self

    def with_bed_file(self, bed_file_records):
        tmp_bed_file = os.path.join(self.__test_case.work_dir, "test.bed")
        # come up with a temporary file name ending with .bed, tmp
        with BEDWriterContextManager(tmp_bed_file) as bed_writer:
            for record in bed_file_records:
                bed_writer.write_bed_record(record)
        self._config['regions'] = tmp_bed_file
        return self

    def with_turn_on_large_variant_calls(self, toggle_value):
        self._config['turnOnLargeVariantCalls'] = toggle_value
        return self

    def with_verbosity(self, verbosity):
        self._config['verbosity'] = verbosity
        return self

    def with_bam_filenames(self, bam_filenames):
        self.__bam_filenames = bam_filenames
        return self

    def call(self, expected_success=True):
        filestem = os.path.join(self.__test_case.work_dir, "_")
        wecall_input_data_builder = WecallInputDataBuilder(
            self.__test_case.work_dir)

        for chrom in self._sample_bank.keys():
            wecall_input_data_builder.with_sample_bank(
                self._sample_bank[chrom])

        if self.__reference_filename is not None:
            wecall_input_data_builder.with_ref_filename(
                self.__reference_filename)

        if self.__bam_filenames is not None:
            wecall_input_data_builder.with_bam_filenames(self.__bam_filenames)

        wecall_input_data = wecall_input_data_builder.build()

        for file_name, (bam_data, with_read_group) in self.__bam_data.items():
            bam_builder = BAMBuilder(
                os.path.join(
                    self.__test_case.work_dir,
                    file_name),
                with_read_group)
            for sample_name, sequence_bank in bam_data.items():
                ref = sequence_bank.reference
                bam_builder.with_bam_contig_data(
                    ref.chrom, ref.length_minus_deletions(), sample_name, sequence_bank)
            bam_builder.build()
            wecall_input_data.bam_filenames.append(bam_builder.filename)

        wecall_config_builder = WecallConfigBuilder(
            wecall_input_data, filestem)
        for key, value in self._config.items():
            wecall_config_builder.with_configuration(key, value)

        vc_wrapper = VariantCallerWrapper(
            filestem, wecall_config_builder.build())

        for bam_filename in wecall_input_data.bam_filenames:
            for chrom, sample_bank in self._sample_bank.items():
                region = "{}:{}".format(chrom, sample_bank.reference.pos_from)
                log_bam_file(
                    wecall_input_data.reference_filename,
                    bam_filename,
                    region)

        if self._output_vcf_filename is not None:
            vc_wrapper.output_vcf = self._output_vcf_filename

        if self.__log_filename is not None:
            vc_wrapper.log_filename = self.__log_filename

        vc_wrapper.run()

        if expected_success:
            self.__test_case.assertEqual(0, vc_wrapper.return_code, vc_wrapper.stderr.decode())
        else:
            self.__test_case.assertNotEqual(0, vc_wrapper.return_code, vc_wrapper.stderr.decode())

        return SVCExpectation(
            self.__test_case,
            vc_wrapper.stdout,
            vc_wrapper.stderr,
            vc_wrapper.return_code,
            vc_wrapper.output_vcf,
            vc_wrapper.log_filename
        )


class SVCExpectation(object):

    def __init__(
            self,
            test_case,
            stdout,
            stderr,
            return_code,
            output_vcf_path,
            log_filename):
        self.__stdout = stdout.decode()
        self.__stderr = stderr.decode()
        self.__return_code = return_code
        self.__output_vcf_path = output_vcf_path
        self.__log_filename = log_filename
        self.__test_case = test_case

    def with_log(self):
        self.__test_case.assertTrue(
            os.path.isfile(
                self.__log_filename),
            self.__log_filename)
        return SVCLogFileExpectation(self.__test_case, self.__log_filename)

    def with_output_vcf(self):
        self.__test_case.assertTrue(
            os.path.isfile(
                self.__output_vcf_path),
            self.__output_vcf_path)
        return VCFExpectation(self.__test_case, self.__output_vcf_path)

    def missing_output_vcf(self):
        self.__test_case.assertFalse(
            os.path.exists(
                self.__output_vcf_path),
            self.__output_vcf_path)
        return self

    def incorrect_var_ids_error(self, *incorrect_var_ids):
        expected_format = " or ".join(
            filter(None, (", ".join(incorrect_var_ids[0:-1]), incorrect_var_ids[-1])))
        self.__test_case.assertRegexpMatches(self.__stderr, re.escape(
            "FAILED - Could not find filter ID(s): " + expected_format + "\n"))
        return self

    def bed_file_contains_contigs_that_are_not_present_in_the_reference_error(
            self, *bad_regions):
        bad_regions_list = list((str(region) for region in bad_regions))
        self.__test_case.assertRegexpMatches(
            self.__stderr,
            re.escape(
                "FAILED - Region(s) " +
                ",".join(bad_regions_list) +
                " are not contained in reference.\n"))
        return self

    def regions_contains_both_bedfile_and_region_string_error(self):
        self.__test_case.assertRegexpMatches(self.__stderr, re.escape(
            "FAILED - Can not have mixture of BED files and region strings\n"))
        return self

    def bedfile_does_not_exist_error(self, bed_file_name):
        self.__test_case.assertRegexpMatches(self.__stderr, re.escape(
            "FAILED - BED file {} does not exist\n".format(bed_file_name)))
        return self

    def genotyping_is_incompatible_with_outputting_reference_calls_error(self):
        self.__test_case.assertRegexpMatches(self.__stderr, re.escape(
            "FAILED - Genotyping is incompatible with outputting reference calls.\n"))
        return self

    def attempt_to_load_invalid_contig_warning(self, contig_name):
        self.__test_case.assertIn(
            "WARNING -- Attempted to load an invalid contig \"{}\" from the BAM file - "
            "Check that the contig names in the reference file match those in the BAM.".format(contig_name),
            self.__stderr)
        return self

    def work_dir_not_a_directory_error(self, location):
        self.__test_case.assertRegexpMatches(
            self.__stderr,
            "FAILED - Working dir: {} is not a directory\n".format(location)
        )
        return self

    def missing_genotype_file(self, filename):
        self.__test_case.assertRegexpMatches(self.__stderr, re.escape(
            "FAILED - Genotype file {} does not exist\n".format(filename)))
        return self

    def output_exists_error(self, filename):
        self.__test_case.assertRegexpMatches(
            self.__stderr, "FAILED - {} already exists".format(filename))
        return self

    def missing_genotype_index_file(self, filename):
        self.__test_case.assertRegexpMatches(self.__stderr, re.escape(
            "FAILED - Genotype index file {} does not exist\n".format(filename)))
        return self

    def unexpected_genotype_file_format(self, filename):
        self.__test_case.assertRegexpMatches(self.__stderr, re.escape(
            "FAILED - File {} does not have .gz extension\n".format(filename)))
        return self

    def with_mem_limit_range_error(self):
        self.__test_case.assertRegexpMatches(
            self.__stderr,
            re.escape("FAILED - <memLimit> not in acceptable range. \n")
        )
        return self

    def with_incorrect_output_format_error(self):
        self.__test_case.assertRegexpMatches(
            self.__stderr, "FAILED - output file format must be VCF4.1 or VCF4.2")
        return self


class SVCLogFileExpectation(object):
    def __init__(self, test_case, path):
        self.__test_case = test_case
        self.__path = path

        with open(self.__path) as fp:
            for line in fp:
                logging.info('> ' + line.rstrip())
            logging.info("")

        self.__lines = open(self.__path).read()

    def bed_file_contains_contigs_that_are_not_present_in_the_reference_warning(
            self, *bad_regions):
        bad_regions_list = list((str(region) for region in bad_regions))
        self.__test_case.assertRegexpMatches(
            self.__lines,
            re.escape(
                "WARNING -- Contig(s) " +
                ",".join(bad_regions_list) +
                " are not contained in reference.\n"))
        return self

    def input_variant_trimmed_warning(self, input, trimmed):
        def we_call_variant(var):
            return "Variant({}:{}-{} {} --> {})".format(var.chrom,
                                                        var.pos_from, var.pos_to, var.ref, var.alt)
        self.__test_case.assertRegexpMatches(
            self.__lines,
            re.escape(
                "WARNING -- Trimming input variant from: " +
                we_call_variant(input) +
                " to " +
                we_call_variant(trimmed) +
                ".\n"))
        return self
