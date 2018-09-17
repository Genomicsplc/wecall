# All content Copyright (C) 2018 Genomics plc
from wecall.bamutils.sequence_bank import SequenceBank
from wecall.genomics.reference_chromosome import ReferenceChromosome
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestBAMWithoutDataOnChromosome(BaseTest):
    def test_should_be_able_to_process_BAM_files_with_missing_data_on_chromosomes(self):
        driver = SVCDriver(self) \
            .with_ref_sequence(
                'CGGCGGTCGAACGGAGCCCCAAGCGAAGCTCAAAACATGG', 0, chrom="1") \
            .with_ref_sequence(
                'CGGCGGTAAAACGGAGCCCCAAGCTTTTTTCAAAACATGG', 0, chrom="2") \
            .with_read(
                '..................T.....................', n_fwd=10, n_rev=10, chrom='1')

        expect = driver.call(expected_success=True)
        expect.attempt_to_load_invalid_contig_warning("2")
        vcf_expect = expect.with_output_vcf()
        vcf_expect.record_count(1)


class TestSingleSampleBAMWithoutReadGroup(BaseTest):
    def test_should_use_sample_name_if_available(self):
        chrom = '14'

        sequence_bank = SequenceBank(ReferenceChromosome(
            'CGGCGGTCGAACGGAGCCCCAAGCGAAGCTCAAAACATGG', 0, chrom))
        sequence_bank.add_sequence(
            '      ...........A.............         ', n_fwd=10, n_rev=10)

        driver = SVCDriver(self).with_ref_sequence(
            'CGGCGGTCGAACGGAGCCCCAAGCGAAGCTCAAAACATGG', chrom=chrom)\
            .with_bam_data('pi.bam', {'sample': sequence_bank}, True)

        expect = driver.call()

        expect.with_output_vcf().record_count(1).with_samples(['sample'])

    def test_should_use_filename_when_no_sample_name_available(self):
        chrom = '14'

        sequence_bank = SequenceBank(ReferenceChromosome(
            'CGGCGGTCGAACGGAGCCCCAAGCGAAGCTCAAAACATGG', 0, chrom))
        sequence_bank.add_sequence(
            '      ...........A.............         ', n_fwd=10, n_rev=10)

        driver = SVCDriver(self).with_ref_sequence(
            'CGGCGGTCGAACGGAGCCCCAAGCGAAGCTCAAAACATGG', chrom=chrom
        ).with_bam_data('pi.bam', {'sample': sequence_bank}, False)

        expect = driver.call()

        expect.with_output_vcf().record_count(1).with_samples(['pi'])
