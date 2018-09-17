# All content Copyright (C) 2018 Genomics plc
import os
import unittest

import wecall.common.exceptions
from wecall.genomics import variant
from wecall.genomics.variant import Variant
from wecall.vcfutils import record
from wecall.vcfutils.genotype_call import GenotypeCall
from wecall.vcfutils.info_data import InfoData
from wecall.vcfutils.parser import VCFReaderContextManager
from wecall.vcfutils.record import common_prefix_length, trimmed_vcf_ref_alt, vcf_row_from_record, \
    variant_from_vcf, variant_quality_from_vcf, variant_ids_from_vcf, \
    vcf_id_entry_from_variant_ids, filters_from_vcf, split_MNP_variant, Record
from wecall.vcfutils.sample_data import SampleData
from wecall.vcfutils.schema import Schema


class TestUtilityFunctionsInRecordModule(unittest.TestCase):
    def test_quality_from_vcf(self):
        quality = "50"
        self.assertEqual(variant_quality_from_vcf(quality), 50.0)
        self.assertEqual(variant_quality_from_vcf("."), None)
        self.assertRaises(wecall.common.exceptions.EchidnaException, variant_quality_from_vcf, "String")

    def test_variant_ids_from_vcf(self):
        self.assertEqual(variant_ids_from_vcf("A,B"), {"A", "B"})
        self.assertEqual(variant_ids_from_vcf("."), set())

    def test_vcf_id_entry_from_variants_ids(self):
        self.assertEqual(vcf_id_entry_from_variant_ids({"A", "B"}), "A,B")
        self.assertEqual(vcf_id_entry_from_variant_ids(set()), ".")

    def test_filters_from_vcf(self):
        self.assertEqual(filters_from_vcf("PASS"), set())
        self.assertEqual(filters_from_vcf("PISS"), {"PISS"})
        self.assertEqual(filters_from_vcf("PASS;PISS"), {"PASS", "PISS"})

    def test_variant_from_vcf(self):
        chrom = "blah"
        pos = 21
        ref = "hell"
        alt = "heaven"
        self.assertEqual(variant_from_vcf(chrom, pos, ref, alt), variant.Variant(chrom, 20, ref, alt))

    def test_split_MNP_variant(self):
        chrom = "blah"
        pos = 21
        ref = "hell"
        alt = "beli"
        var = variant.Variant(chrom, pos, ref, alt)
        self.assertEqual(
            list(split_MNP_variant(var)),
            [
                variant.Variant(chrom, pos, ref[0], alt[0]),
                variant.Variant(chrom, pos + 3, ref[3], alt[3]),
            ]
        )
        self.assertEqual(
            list(split_MNP_variant(var, include_ref_calls=True)),
            [
                variant.Variant(chrom, pos, ref[0], alt[0]),
                variant.Variant(chrom, pos + 1, ref[1], alt[1]),
                variant.Variant(chrom, pos + 2, ref[2], alt[2]),
                variant.Variant(chrom, pos + 3, ref[3], alt[3]),
            ]
        )


class RecordTest(unittest.TestCase):
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "example_data")
        self.work_dir = os.path.join(os.environ["ECHIDNA_TEST_RESULTS"], *self.id().split("."))
        try:
            os.makedirs(self.work_dir)
        except OSError:
            pass

    def __get_example_schema(self, filename):
        with VCFReaderContextManager(os.path.join(self.data_dir, filename)) as vcf_handler:
            vcf_handler.read_header()
            return vcf_handler.header

    def test_eq(self):
        reference = Record(None, Variant("1", 20, "A", "G"), set(), 0.0, set(),
                           InfoData(None, {}), SampleData([], []), False)

        self.assertTrue(reference == Record(None, Variant("1", 20, "A", "G"), set(), 0.0, set(),
                                            InfoData(None, {}), SampleData([], []), False))

        self.assertFalse(reference == Record(None, Variant("2", 20, "A", "G"), set(), 0.0, set(),
                                             InfoData(None, {}), SampleData([], []), False))

        self.assertFalse(reference == Record(None, Variant("1", 20, "A", "G"), set("rs0"), 0.0, set(),
                                             InfoData(None, {}), SampleData([], []), False))

        self.assertFalse(reference == Record(None, Variant("1", 20, "A", "G"), set(), 5.0, set(),
                                             InfoData(None, {}), SampleData([], []), False))

        self.assertFalse(reference == Record(None, Variant("1", 20, "A", "G"), set(), 0.0, set("CV"),
                                             InfoData(None, {}), SampleData([], []), False))

        self.assertFalse(reference == Record(None, Variant("1", 20, "A", "G"), set(), 0.0, set(),
                                             InfoData(None, {'AF': []}), SampleData([], []), False))

        self.assertFalse(reference == Record(None, Variant("1", 20, "A", "G"), set(), 0.0, set(),
                                             InfoData(None, {}), SampleData([], ['NA12787']), False))

        self.assertFalse(reference == Record(None, Variant("1", 20, "A", "G"), set(), 0.0, set(),
                                             InfoData(None, {}), SampleData([], []), True))

    def test_ne(self):
        reference = Record(None, Variant("1", 20, "A", "G"), set(), 0.0, set(),
                           InfoData(None, {}), SampleData([], []), False)

        self.assertFalse(reference != Record(None, Variant("1", 20, "A", "G"), set(), 0.0, set(),
                                             InfoData(None, {}), SampleData([], []), False))

        self.assertTrue(reference != Record(None, Variant("2", 20, "A", "G"), set(), 0.0, set(),
                                            InfoData(None, {}), SampleData([], []), False))

        self.assertTrue(reference != Record(None, Variant("1", 20, "A", "G"), set("rs0"), 0.0, set(),
                                            InfoData(None, {}), SampleData([], []), False))

        self.assertTrue(reference != Record(None, Variant("1", 20, "A", "G"), set(), 5.0, set(),
                                            InfoData(None, {}), SampleData([], []), False))

        self.assertTrue(reference != Record(None, Variant("1", 20, "A", "G"), set(), 0.0, set("CV"),
                                            InfoData(None, {}), SampleData([], []), False))

        self.assertTrue(reference != Record(None, Variant("1", 20, "A", "G"), set(), 0.0, set(),
                                            InfoData(None, {'AF': []}), SampleData([], []), False))

        self.assertTrue(reference != Record(None, Variant("1", 20, "A", "G"), set(), 0.0, set(),
                                            InfoData(None, {}), SampleData([], ['NA12787']), False))

        self.assertTrue(reference != Record(None, Variant("1", 20, "A", "G"), set(), 0.0, set(),
                                            InfoData(None, {}), SampleData([], []), True))

    def test_read_variant_from_vcf(self):
        with VCFReaderContextManager(os.path.join(self.data_dir, "vcf_example.vcf")) as vcf_handler:
            variant_gen = (record.variant for record in vcf_handler.read_records())
            next_variant = next(variant_gen)
            self.assertEqual(next_variant.chrom, "20")
            self.assertEqual(next_variant.pos_from, 9)
            self.assertEqual(next_variant.ref, "CT")
            self.assertEqual(next_variant.alt, "C")

    def test_read_record_line(self):
        with VCFReaderContextManager(os.path.join(self.data_dir, "vcf_example.vcf")) as vcf_handler:

            record_gen = vcf_handler.read_records()
            next_record = next(record_gen)

            self.assertEqual(next_record.chrom, "20")
            self.assertEqual(next_record.pos_from, 9)
            self.assertEqual(next_record.ids, set())
            self.assertEqual(next_record.ref, "CT")
            self.assertEqual(next_record.alt, "C")
            self.assertEqual(next_record.quality, 3000)
            self.assertEqual(next_record.filters, set())
            self.assertEqual(next_record.passes_filter, True)
            self.assertEqual(next_record.from_multi_alt, False)
            self.assertEqual(next_record.type, variant.TYPE_DEL)

            self.assertEqual(next_record.info['PP'], [3000])
            self.assertEqual(next_record.info['DP'], [250])
            self.assertEqual(next_record.info['VC'], [100])
            self.assertEqual(next_record.info['ABPV'], [0.2])
            self.assertEqual(next_record.info['SBPV'], [0.3])
            self.assertEqual(next_record.info['MQ'], [70])
            self.assertEqual(next_record.info['QD'], [None])

            self.assertTrue(next_record.sample_info.has_sample("sample1"))
            self.assertEqual(next_record.genotypes, {"sample1": GenotypeCall("1|0"), "sample2": GenotypeCall("1|1")})
            self.assertEqual(next_record.sample_info.get_field("sample1", 'GT'), GenotypeCall("1|0"))
            self.assertEqual(next_record.sample_info.get_field("sample1", 'PL'), [3000, 0, 3000])
            self.assertEqual(next_record.sample_info.get_field("sample1", "GQ"), [1000])

    def test_read_sample_data(self):
        schema = self.__get_example_schema("vcf_example.vcf")
        sample_schema = [key for key, _ in schema.iter_sample_data()]

        sample_data = SampleData(sample_schema, ['sample1'])

        sample_data.add_sample_data("sample1", "GT", GenotypeCall("1|0"))
        sample_data.add_sample_data("sample1", "PL", [3000, 0, 3000])
        sample_data.add_sample_data("sample1", "GQ", [1000])
        sample_data.add_sample_data("sample1", "PQ", [2000])
        sample_data.add_sample_data("sample1", "PS", [60000])
        sample_data.add_sample_data("sample1", "AD", [140, 110])
        sample_data.add_sample_data("sample1", "DP", [250])
        sample_data.add_sample_data("sample1", "VAF", [0.4])

        self.assertTrue(sample_data.has_sample("sample1"))
        self.assertEqual(sample_data.genotypes(), {"sample1": GenotypeCall("1|0")})
        self.assertEqual(sample_data.get_field("sample1", 'GT'), GenotypeCall("1|0"))
        self.assertEqual(sample_data.get_field("sample1", 'PL'), [3000, 0, 3000])

        genotype_data = sample_data.get_genotype_data("sample1")
        self.assertEqual(genotype_data.genotype(), GenotypeCall("1|0"))
        self.assertEqual(genotype_data['GT'], GenotypeCall("1|0"))
        self.assertEqual(genotype_data['PL'], [3000, 0, 3000])

    def test_should_write_missing_values_in_sample_data(self):
        with VCFReaderContextManager(os.path.join(self.data_dir, "vcf_example.vcf")) as vcf_handler:
            first_record = next(vcf_handler.read_records())

        sample_data = SampleData(['GT', 'PL', 'GQ'], ['sample1', 'sample2', 'sample3'])

        sample_data.add_sample_data("sample1", "GT", GenotypeCall("1|0"))
        sample_data.add_sample_data("sample1", "PL", [3000, 0, 3000])
        sample_data.add_sample_data("sample1", "GQ", [1000])

        sample_data.add_sample_data("sample2", "GT", GenotypeCall("1|1"))
        sample_data.add_sample_data("sample2", "PL", [2000, 0, 1000])
        sample_data.add_sample_data("sample2", "GQ", [3])

        first_record.sample_info = sample_data

        print((sample_data.to_vcf_columns()))
        vcf_string = vcf_row_from_record(first_record)
        expected_vcf_string = "20	10	.	CT	C	3000	PASS	PP=3000;DP=250;DPR=140;DPF=110;VC=100;VCR=49;VCF=51;ABPV=0.2;SBPV=0.3;MQ=70.0;BR=31.0;QD=None	GT:PL:GQ	1|0:3000,0,3000:1000	1|1:2000,0,1000:3	./.:.:."  # noqa
        self.assertEqual(expected_vcf_string, vcf_string)

    def test_should_return_default_diploid_genotype(self):
        sample_data = SampleData(['GT', 'GL'], ["NA12878"])

        self.assertEqual(GenotypeCall("./."), GenotypeCall("./."))

        self.assertTrue(sample_data.has_sample("NA12878"))
        self.assertEqual(sample_data.genotypes(), {"NA12878": GenotypeCall("./.")})
        self.assertEqual(sample_data.get_field("NA12878", 'GT'), GenotypeCall("./."))
        self.assertEqual(sample_data.get_field("NA12878", 'GL'), [])

        genotype_data = sample_data.get_genotype_data("NA12878")
        self.assertEqual(genotype_data.genotype(), GenotypeCall("./."))
        self.assertEqual(genotype_data['GT'], GenotypeCall("./."))
        self.assertEqual(genotype_data['GL'], [])

    def test_split_empty_sample_data_string(self):
        schema = self.__get_example_schema("vcf_example.vcf")
        cols = """1\t11082325\tRS1\tG\tC,A\t.\t.\tPP=.;DP=.;DPR=.;DPF=.;VC=.;VCR=.;VCF=.;ABPV=.;SBPV=.;MQ=.;BR=.;QD=.\tGT:PL:GQ\t1|0:3000,0,3000:1000\t1|1:2000,0,1000:3""".split()  # noqa
        first_record = next(record.generate_records(schema, cols))
        self.assertEqual(first_record.alt, 'C')
        self.assertTrue(first_record.from_multi_alt)

        second_record = next(record.generate_records(schema, cols))
        self.assertEqual(first_record.info, second_record.info)

    def test_split_genotype_calls(self):
        self.assertEqual(record.split_GT("0/0", 2), ["0/0", "0/0"])
        self.assertEqual(record.split_GT("0/1", 2), ["0/1", "0/0"])
        self.assertEqual(record.split_GT("0/2", 2), ["0/0", "0/1"])
        self.assertEqual(record.split_GT("1/1", 2), ["1/1", "0/0"])
        self.assertEqual(record.split_GT("1/2", 2), ["0/1", "0/1"])
        self.assertEqual(record.split_GT("2/2", 2), ["0/0", "1/1"])

        self.assertEqual(record.split_GT("0/0", 3), ["0/0", "0/0", "0/0"])
        self.assertEqual(record.split_GT("0/1", 3), ["0/1", "0/0", "0/0"])
        self.assertEqual(record.split_GT("0/2", 3), ["0/0", "0/1", "0/0"])
        self.assertEqual(record.split_GT("0/3", 3), ["0/0", "0/0", "0/1"])
        self.assertEqual(record.split_GT("1/1", 3), ["1/1", "0/0", "0/0"])
        self.assertEqual(record.split_GT("1/2", 3), ["0/1", "0/1", "0/0"])
        self.assertEqual(record.split_GT("1/3", 3), ["0/1", "0/0", "0/1"])
        self.assertEqual(record.split_GT("2/2", 3), ["0/0", "1/1", "0/0"])
        self.assertEqual(record.split_GT("2/3", 3), ["0/0", "0/1", "0/1"])
        self.assertEqual(record.split_GT("3/3", 3), ["0/0", "0/0", "1/1"])

    def test_split_unknown_genotype_calls(self):
        self.assertEqual(record.split_GT("./.", 2), ["./.", "./."])
        self.assertEqual(record.split_GT("./1", 2), ["./1", "./0"])
        self.assertEqual(record.split_GT("./2", 2), ["./0", "./1"])

    def test_split_phased_genotype_calls(self):
        self.assertEqual(record.split_GT("0|0", 2), ["0|0", "0|0"])
        self.assertEqual(record.split_GT("0|1", 2), ["0|1", "0|0"])
        self.assertEqual(record.split_GT("0|2", 2), ["0|0", "0|1"])
        self.assertEqual(record.split_GT("1|0", 2), ["1|0", "0|0"])
        self.assertEqual(record.split_GT("1|1", 2), ["1|1", "0|0"])
        self.assertEqual(record.split_GT("1|2", 2), ["1|0", "0|1"])
        self.assertEqual(record.split_GT("2|0", 2), ["0|0", "1|0"])
        self.assertEqual(record.split_GT("2|1", 2), ["0|1", "1|0"])
        self.assertEqual(record.split_GT("2|2", 2), ["0|0", "1|1"])

        self.assertEqual(record.split_GT("0|0", 3), ["0|0", "0|0", "0|0"])
        self.assertEqual(record.split_GT("0|1", 3), ["0|1", "0|0", "0|0"])
        self.assertEqual(record.split_GT("0|2", 3), ["0|0", "0|1", "0|0"])
        self.assertEqual(record.split_GT("0|3", 3), ["0|0", "0|0", "0|1"])
        self.assertEqual(record.split_GT("1|0", 3), ["1|0", "0|0", "0|0"])
        self.assertEqual(record.split_GT("1|1", 3), ["1|1", "0|0", "0|0"])
        self.assertEqual(record.split_GT("1|2", 3), ["1|0", "0|1", "0|0"])
        self.assertEqual(record.split_GT("1|3", 3), ["1|0", "0|0", "0|1"])
        self.assertEqual(record.split_GT("2|0", 3), ["0|0", "1|0", "0|0"])
        self.assertEqual(record.split_GT("2|1", 3), ["0|1", "1|0", "0|0"])
        self.assertEqual(record.split_GT("2|2", 3), ["0|0", "1|1", "0|0"])
        self.assertEqual(record.split_GT("2|3", 3), ["0|0", "1|0", "0|1"])
        self.assertEqual(record.split_GT("3|0", 3), ["0|0", "0|0", "1|0"])
        self.assertEqual(record.split_GT("3|1", 3), ["0|1", "0|0", "1|0"])
        self.assertEqual(record.split_GT("3|2", 3), ["0|0", "0|1", "1|0"])
        self.assertEqual(record.split_GT("3|3", 3), ["0|0", "0|0", "1|1"])

    def __cleanup_tmp_files(self, file_name):
        try:
            os.remove(file_name)
        except OSError as ex:
            self.assertTrue(False, str(ex))
        else:
            pass


class TestCommonPrefixLength(unittest.TestCase):

    def test_no_common_prefix(self):
        self.assertEqual(0, common_prefix_length("ACDE", "BCDE"))

    def test_common_prefix_matching_strings(self):
        self.assertEqual(4, common_prefix_length("ABCA", "ABCA"))

    def test_common_prefix_then_end_of_one_string(self):
        self.assertEqual(2, common_prefix_length("AAA", "AA"))
        self.assertEqual(2, common_prefix_length("AA", "AAA"))

    def test_common_prefix_different_suffix(self):
        self.assertEqual(3, common_prefix_length("AAAA", "AAAB"))
        self.assertEqual(3, common_prefix_length("AAAB", "AAAA"))


class TestTrimmedVCFRefAlt(unittest.TestCase):

    def test_refcall(self):
        self.assertEqual(trimmed_vcf_ref_alt("A", "A"), (0, "A", "A"))

    def test_fails_with_long_refcall(self):
        self.assertRaises(wecall.common.exceptions.EchidnaException, trimmed_vcf_ref_alt, "ACDE", "ACDE")

    def test_with_SNP(self):
        self.assertEqual(trimmed_vcf_ref_alt("ABCDE", "ABQDE"), (2, "C", "Q"))

    def test_with_trimmed_SNP(self):
        self.assertEqual(trimmed_vcf_ref_alt("C", "Q"), (0, "C", "Q"))

    def test_with_MNP(self):
        self.assertEqual(trimmed_vcf_ref_alt("ABCDEFG", "ABQDPFG"), (2, "CDE", "QDP"))

    def test_with_trimmed_MNP(self):
        self.assertEqual(trimmed_vcf_ref_alt("ABCDEFG", "XBCDEFY"), (0, "ABCDEFG", "XBCDEFY"))

    def test_with_INS(self):
        self.assertEqual(trimmed_vcf_ref_alt("ABCDEF", "ABCPQDEF"), (2, "C", "CPQ"))

    def test_with_INS_at_start(self):
        self.assertEqual(trimmed_vcf_ref_alt("DEF", "PQDEF"), (0, "D", "PQD"))

    def test_with_DEL(self):
        self.assertEqual(trimmed_vcf_ref_alt("ABCPQDEF", "ABCDEF"), (2, "CPQ", "C"))

    def test_with_DEL_at_start(self):
        self.assertEqual(trimmed_vcf_ref_alt("PQDEF", "DEF"), (0, "PQD", "D"))

    def test_with_common_prefix_and_suffix(self):
        self.assertEqual(trimmed_vcf_ref_alt("AT", "ATAT"), (0, "A", "ATA"))

    def test_fails_with_monomorphic_variant(self):
        self.assertRaises(wecall.common.exceptions.EchidnaException, trimmed_vcf_ref_alt, "AT", ".")
        self.assertRaises(wecall.common.exceptions.EchidnaException, trimmed_vcf_ref_alt, ".", "AT")

    def test_fails_with_empty_ref(self):
        self.assertRaises(wecall.common.exceptions.EchidnaException, trimmed_vcf_ref_alt, "", "A")

    def test_fails_with_empty_alt(self):
        self.assertRaises(wecall.common.exceptions.EchidnaException, trimmed_vcf_ref_alt, "A", "")

    def test_fails_with_empty_ref_and_alt(self):
        self.assertRaises(wecall.common.exceptions.EchidnaException, trimmed_vcf_ref_alt, "", "")


class TestInfoFieldFormatting(unittest.TestCase):
    def test_should_format_a_present_flag(self):
        schema = Schema()
        schema.set_info_data('F', '0', 'Flag', 'Flag')
        info_data = InfoData(schema, {"F": None})
        self.assertEqual('F', info_data.to_vcf())

    def test_should_format_no_data(self):
        info_data = InfoData(None, {})
        self.assertEqual('.', info_data.to_vcf())

    def test_should_format_a_string(self):
        info_data = InfoData(None, {'K': 'V'})
        self.assertEqual('K=V', info_data.to_vcf())

    def test_should_format_a_string_list(self):
        schema = Schema()
        schema.set_info_data('K', 'A', 'String', 'K')
        info_data = InfoData(schema, {'K': ['V1', 'V2']})
        self.assertEqual('K=V1,V2', info_data.to_vcf())

    def test_should_format_an_int_list(self):
        schema = Schema()
        schema.set_info_data('K', 'A', 'Integer', 'K')
        info_data = InfoData(schema, {'K': [1, 2, 3]})
        self.assertEqual('K=1,2,3', info_data.to_vcf())

    def test_should_format_a_float_list(self):
        schema = Schema()
        schema.set_info_data('K', 'A', 'Integer', 'K')
        info_data = InfoData(schema, {'K': [1.0, 2.66, 3.0]})
        self.assertEqual('K=1.0,2.66,3.0', info_data.to_vcf())

    def test_should_format_multiple_values(self):
        schema = Schema()
        schema.set_info_data('K', 'A', 'Float', 'K')
        schema.set_info_data('K2', 'A', 'String', 'K')
        schema.set_info_data('K3', '0', 'Flag', 'K')
        schema.set_info_data('K4', 'A', 'String', 'K')
        info_data = InfoData(schema, {'K3': None, 'K2': ['S2'], 'K': [1.0, 2.66, 3.0], 'K4': ['S4']})
        self.assertEqual('K=1.0,2.66,3.0;K2=S2;K3;K4=S4', info_data.to_vcf())
