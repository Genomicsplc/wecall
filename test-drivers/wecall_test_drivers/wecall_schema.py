# All content Copyright (C) 2018 Genomics plc
from wecall.vcfutils.schema import Schema


def wecall_schema(
        file_date=None,
        reference=None,
        contigs=None,
        add_ref_calls=True,
        format='4.2'):
    schema = Schema()
    if file_date is not None:
        schema.file_metadata['fileDate'] = file_date
    if reference is not None:
        schema.file_metadata['reference'] = reference

    app_name = 'weCall'
    version_number = '2.0.1'
    app = {'4.1': None, '4.2': app_name}[format]
    version = {'4.1': None, '4.2': version_number}[format]

    schema.file_metadata[
        'disclaimer'] = 'This software is in beta-testing. Results generated using the software are confidential and should only be used for research purposes in accordance with the legal agreement with Genomics plc.'  # noqa
    schema.file_metadata['source'] = '{application!s} v{version!s}'.format(
        application=app_name, version=version_number)  # noqa

    schema.set_info_data(
    'ABPV',
    'A',
    'Float',
    'Allele bias P-value; probability that fraction of reads supporting alt allele (VC) amongst read depth (DP) is '
    'more extreme than expected assuming a beta-binomial distribution.',
    app,
     version)  # noqa
    schema.set_info_data(
    'MQ',
    'A',
    'Float',
    'Root mean square of mapping quality of reads supporting each alternative allele.',
    app,
     version)  # noqa
    schema.set_info_data(
    'PP',
    'A',
    'Integer',
    'Posterior probability (phred scaled) that this variant does not segregate.',
    app,
     version)  # noqa
    schema.set_info_data(
    'SBPV',
    'A',
    'Float',
    'Strand bias P-value; probability that the fraction of forward reads (VCF) amongst reads supporting alt allele '
    '(VC) is more extreme than expected assuming a beta-binomial distribution.',
    app,
     version)  # noqa
    schema.set_info_data(
        'DP',
        '1',
        'Integer',
        'Total depth of read coverage at this locus.',
        app,
        version)
    schema.set_info_data(
    'DPF',
    '1',
    'Integer',
    'Total probabilistic depth of forward read coverage at this locus (sum of probabilities of each read supporting '
    'the variant).',
    app,
     version)  # noqa
    schema.set_info_data(
    'DPR',
    '1',
    'Integer',
    'Total probabilistic depth of reverse read coverage at this locus (sum of probabilities of each read supporting '
    'the variant).',
    app,
     version)  # noqa
    schema.set_info_data(
    'VC',
    'A',
    'Integer',
    'Total probabilistic number of reads supporting each alternative allele (sum of probabilities of each read '
    'supporting the allele).',
    app,
     version)  # noqa
    schema.set_info_data(
    'VCF',
    'A',
    'Integer',
    'Total probabilistic number of forward reads supporting each alternative allele (sum of probabilities of '
    'each read supporting the allele).',
    app,
     version)  # noqa
    schema.set_info_data(
    'VCR',
    'A',
    'Integer',
    'Total probabilistic number of reverse reads supporting each alternative allele (sum of probabilities of each '
    'read supporting the allele).',
    app,
     version)  # noqa
    schema.set_info_data(
    'QD',
    'A',
    'Float',
    'Ratio of phred-scaled posterior probability (PP) to number of supporting reads for each allele (VC).',
    app,
     version)  # noqa
    schema.set_info_data(
    'BR',
    'A',
    'Float',
    'The median of the per-read min base quality (within a interval of the locus) taken over reads supporting '
    'each allele.',
    app,
     version)  # noqa

    schema.set_sample_data(
        'GT',
        '1',
        'String',
        'Genotypes of reference and alternative alleles in order listed.')

    if add_ref_calls:
        schema.set_info_data(
            'BEG',
            '1',
            'Integer',
            'Start position of reference call block.',
            app,
            version)
        schema.set_info_data(
            'END',
            '1',
            'Integer',
            'End position of reference call block (inclusive).',
            app,
            version)
        schema.set_info_data(
            'LEN',
            '1',
            'Integer',
            'Length of reference call block.',
            app,
            version)
        schema.set_sample_data(
            'MIN_DP',
            '1',
            'Integer',
            'Minimum read coverage observed within the reference block.')

    schema.set_sample_data('GQ', '1', 'Integer',
                           'Phred-scaled genotype quality (i.e. posterior probability that the genotype call is incorrect).')  # noqa
    schema.set_sample_data('PQ', '1', 'Integer',
                           'Phred-scaled phase quality (i.e. posterior probability that the phasing is incorrect).')  # noqa
    schema.set_sample_data('PS', '1', 'String', 'Phase set id.')  # noqa
    schema.set_sample_data('PL', 'G', 'Integer',
                           "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification.")  # noqa
    schema.set_sample_data('DP', '1', 'Integer',
                           'Number of reads overlapping the variant site (i.e. INFO::DP split out by sample). For reference calls the average depth (rounded to the nearest integer) over the region is reported.')  # noqa
    schema.set_sample_data('AD', '.', 'Integer',
                           'Probabilistic allelic depths for the ref and alt alleles in the order listed (i.e. INFO::VC split out by sample).')  # noqa
    schema.set_sample_data('VAF', 'A', 'Float',
                           'Probabilistic variant allelic frequencies for each alt allele (FORMAT::AD / FORMAT::DP).')  # noqa

    schema.set_filter('AB',
                      'Allele Bias: Indicates lower number of reads supporting variant than expected (any of INFO::ABPV < 0.009).')  # noqa
    schema.set_filter('SB',
                      'Strand Bias: Indicates imbalance between number of forward and reverse reads supporting variant (any of INFO::SBPV < 0.01).')  # noqa
    schema.set_filter('AB+SB',
                      'Allele + Strand Bias: Indicates that both the AB and SB filters are close to being triggered (any of INFO::ABPV + INFO::SBPV < 0.07).')  # noqa
    schema.set_filter('MQ',
                      'low Mapping Quality: Indicates presence of low mapping quality (any of INFO::MQ < 25).')  # noqa
    schema.set_filter('QD',
                      'Quality over Depth: Indicates low quality relative to number of supporting reads (any of INFO::QD < 3.5 for Indels or INFO::QD < 8 otherwise).')  # noqa
    schema.set_filter('BR',
                      'Bad Reads: Indicates low quality base pairs on reads in the vicinity of variant locus (any of INFO::BR < 0).')  # noqa
    schema.set_filter('NC', 'Not called: Indicates a variant that was not positively genotyped in any sample.')  # noqa
    schema.set_filter('LQ', 'Low Quality: Indicates a low variant quality (any of INFO::PP < 10).')  # noqa

    if contigs is not None:
        for contig_name, contig_data in contigs.items():
            schema.set_contig(contig_name, **contig_data)
    return schema
