# All content Copyright (C) 2018 Genomics plc
import wecall.genomics.variant
import wecall.vcfutils.parser


class VariantCallSet(object):
    """
    Variant Call container to be used ONLY when call set is small. Not in benchmark.
    """

    def __init__(self, test_case):
        self.__test_case = test_case
        self.__vcf_header = None
        self.__vcf_variants = {}

    def __len__(self):
        return len(self.__vcf_variants)

    def get_variants(self):
        return set(self.__vcf_variants.keys())

    def get_variant_records(self):
        return self.__vcf_variants

    def get_variants_with_genotypes(self):
        variants_with_genotypes = {}
        for variant, record in self.__vcf_variants.items():
            variants_with_genotypes[variant] = record.genotypes

        return variants_with_genotypes

    def add_vcf_variants(self, output_vcf_path):
        with wecall.vcfutils.parser.VCFReaderContextManager(output_vcf_path) as vcf_reader:
            vcf_reader.read_header()
            self.__vcf_header = vcf_reader.header
            for record in vcf_reader.read_records():
                self.__test_case.assertFalse(
                    record.variant in self.__vcf_variants,
                    "Repeated variants in VCF.")
                self.__vcf_variants[record.variant] = record

        # ensure genotype likelihoods are within range throughout
        for record in self.__vcf_variants.values():
            for sample_name in record.sample_info.get_sample_names():
                try:
                    for GL_value in record.sample_info.get_field(
                            sample_name, 'GL'):
                        self.__test_case.assertTrue(GL_value <= 0.0)
                except KeyError:
                    pass

    def expect_info_header(
            self,
            ID,
            number=None,
            data_type=None,
            description=None):
        self.__test_case.assertIn(
            ID, {key for key, _ in self.__vcf_header.iter_info_data()})
        item = self.__vcf_header.get_info_data(ID)
        if number is not None:
            self.__test_case.assertEqual(number, item.number)
        if data_type is not None:
            self.__test_case.assertEqual(data_type, item.data_type)
        if description is not None:
            self.__test_case.assertEqual(description, item.description)
        return self

    def assertVariantNotCalled(self, chrom, vcf_pos, ref, alt):
        expected_variant = self.__build_variant(chrom, vcf_pos, ref, alt)
        self.__test_case.assertNotIn(expected_variant, self.__vcf_variants)

    def __assertVariantCalled(self, chrom, vcf_pos, ref, alt):
        expected_variant = self.__build_variant(chrom, vcf_pos, ref, alt)
        self.__test_case.assertIn(expected_variant, self.__vcf_variants)
        return self.__vcf_variants[expected_variant]

    def expect_variant(
            self,
            variant,
            filters=None,
            variant_ids=None,
            require_info_annotations=None,
            require_info_annotation_keys=None,
            missing_info_annotations=None,
    ):
        self.assertVariantCalledWithMetadata(
            variant.chrom,
            variant.one_indexed_pos_from,
            variant.ref,
            variant.alt,
            filters,
            variant_ids,
            require_info_annotations,
            require_info_annotation_keys,
            missing_info_annotations)
        return self

    def assertVariantCalledWithMetadata(
            self, chrom, vcf_pos, ref, alt,
            filters=None,
            variant_ids=None,
            require_info_annotations=None,
            require_info_annotation_keys=None,
            missing_info_annotations=None,
    ):
        generated_record = self.__assertVariantCalled(chrom, vcf_pos, ref, alt)
        if filters is not None:
            self.__check_filters(generated_record, filters)
        if variant_ids is not None:
            self.__check_id(generated_record, variant_ids)
        if require_info_annotation_keys is not None:
            self.__check_required_annotation_keys(
                generated_record.info, require_info_annotation_keys)
        if require_info_annotations is not None:
            self.__check_required_annotations(
                generated_record.info, require_info_annotations)
        if missing_info_annotations is not None:
            self.__check_missing_annotations(
                generated_record.info, missing_info_annotations)

    def assertVariantCalledWithSampleMetadata(
            self, chrom, vcf_pos, ref, alt,
            sample_name,
            require_format_annotations=None,
            require_format_annotation_keys=None,
            missing_format_annotations=None
    ):
        generated_record = self.__assertVariantCalled(chrom, vcf_pos, ref, alt)
        self.__test_case.assertTrue(
            generated_record.sample_info.has_sample(sample_name))
        sample_data = generated_record.sample_info.get_genotype_data(
            sample_name)

        if require_format_annotation_keys is not None:
            self.__check_required_annotation_keys(
                sample_data, require_format_annotation_keys)
        if require_format_annotations is not None:
            self.__check_required_annotations(
                sample_data, require_format_annotations)
        if missing_format_annotations is not None:
            self.__check_missing_annotations(
                sample_data, missing_format_annotations)

    def __check_filters(self, record, filters):
        self.__test_case.assertEqual(record.filters, set(filters))

    def __check_id(self, record, variant_ids):
        self.__test_case.assertEqual(record.ids, variant_ids)
        if variant_ids != set():
            self.__test_case.assertIn("CV", record.info)

    def __check_required_annotations(self, data, annotations):
        self.__test_case.assertDictContainsSubset(annotations, data)

    def __check_missing_annotations(self, data, annotations):
        for item in annotations:
            self.__test_case.assertNotIn(item, data)

    def __check_required_annotation_keys(self, data, annotations):
        for item in annotations:
            self.__test_case.assertIn(item, data)

    def __build_variant(self, chrom, vcf_pos, ref, alt):
        zero_indexed_pos = vcf_pos - 1
        return wecall.genomics.variant.Variant(
            chrom, zero_indexed_pos, ref, alt)
