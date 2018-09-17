# All content Copyright (C) 2018 Genomics plc
import pysam

RG_ID = "Test"


class BAMBuilder(object):

    def __init__(self, filename, with_read_group=True):
        self.__bam_data = None
        self.__filename = filename
        self.__with_read_group = with_read_group

    def with_bam_contig_data(
            self,
            chrom,
            chrom_length,
            sample_name,
            sequence_bank):
        if self.__bam_data is None:
            self.__bam_data = []
        self.__bam_data.append(
            BAMContigData(
                chrom,
                chrom_length,
                sample_name,
                sequence_bank))
        return self

    @property
    def filename(self):
        return self.__filename

    @property
    def header(self):
        sample_names = {
            contig_data.sample_name for contig_data in self.__bam_data}
        bam_header = {'HD': {'VN': '1.0'}, 'SQ': [
            {'LN': contig.chrom_length, 'SN': contig.chrom} for contig in self.__bam_data], }
        if self.__with_read_group:
            bam_header['RG'] = [{"ID": RG_ID + "_" + sample_name,
                                 "SM": sample_name} for sample_name in sorted(sample_names)]
        return bam_header

    def build(self):
        bam_fp = pysam.Samfile(self.__filename, "wb", header=self.header)

        for contig_data in self.__bam_data:
            chrom_id = bam_fp.gettid(contig_data.chrom)
            read_tags = []
            if self.__with_read_group:
                read_tags.append(("RG", RG_ID + "_" + contig_data.sample_name))

            # Sort sequences before writing to bam file. Indexing will
            # otherwise fail.
            reads = sorted(
                contig_data.sequence_bank.build_reads(chrom_id, read_tags),
                key=lambda x: (x.pos, x.seq, x.qual, x.cigarstring, x.mapq)
            )

            for read in reads:
                bam_fp.write(read)

        bam_fp.close()
        pysam.index(self.filename)


class BAMContigData(object):

    def __init__(self, chrom, chrom_length, sample_name, sequence_bank):
        self.chrom = chrom
        self.chrom_length = chrom_length
        self.sample_name = sample_name
        self.sequence_bank = sequence_bank
