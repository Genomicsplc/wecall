# All content Copyright (C) 2018 Genomics plc
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest
from wecall_test_drivers.svc_driver import SVCDriver


ref_alt = "<NON_REF>"


class TestRefCallingMaxRefCallSize(AsciiWecallRunnerTest):
    def test_splits_reference_call_into_three_records(self):
        chrom = "1"
        sample = "bah.asdhaslkdghalsdkfq25451c`52980biqweuo8!"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", chrom=chrom
        ).with_read(
            "                                         ", sample_name=sample
        ).with_output_ref_calls(True).with_max_ref_call_size(20)

        vcf_expect = driver.call().with_output_vcf()

        vcf_expect.record_count(3)
        vcf_expect.has_record(chrom, 0, "A", ref_alt).with_sample(sample).has_genotype("0/0")
        vcf_expect.has_record(chrom, 20, "A", ref_alt).with_sample(sample).has_genotype("0/0")
        vcf_expect.has_record(chrom, 40, "A", ref_alt).with_sample(sample).has_genotype("0/0")
