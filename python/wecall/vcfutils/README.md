# VCF File Parser

# Usage

## VCFReader

Create a VCFReader and read `records`.
`record` is a representation of a single variant call in the vcf file.  Lines with multiple alts are split into single alt `records`.  The method `read_records()` provides a record generator.

Short example (read all variant calls and print position of all SNPs to the screen)

```python
from echidna.vcfutils.parser import VCFReaderContextManager
from echidna.genomics import variant

with VCFReaderContextManager('myFile.vcf') as vcf_reader:
    for record in vcf_reader.read_records():
        if record.type == variant.TYPE_SNP:
            print(record.pos_from)
```

### Available `record` and `variant` properties:

* `chrom` - chrom string (also available in `variant`)
* `pos_from` - integer zero-indexed position (for internal use) (also available in `variant`)
* `one_indexed_pos_from` - integer one-indexed position as used in the wide world (and in vcf) (also available in `variant`)
* `ref` -  (also available in `variant`)
* `alt` -  (also available in `variant`)
* `type` - available types: TYPE_SNP , TYPE_INS, TYPE_DEL, TYPE_MNP -  (also available in `variant`)
* `info` - returns a dictionary for the INFO field
* `sample_info` - returns a dictionary of dictionaries keyed on sample name.  Internal dictionaries are keyed on the properties.  I.e. the 'GL' field for the 'NA12878' can be retrieved by
`record.sample_info['NA12878']['GL']`
* `from_multi_alt` - flag for records that come from a multiallelic line in the .vcf file

For full list of available record properties see `vcfutils/record.py` and `genomics/variant.py`

## Accessing Sample Info Fields:
Sample info fields are accessable as dictionary.
Example: Print all lines that have a homozygous variant for the sample "yoursample".

```python
from echidna.vcfutils.parser import VCFReaderContextManager
with VCFReaderContextManager('inFile.vcf') as vcf_reader:
    for record in vcf_reader.read_records():
        if record.sample_info['yoursample']['GT'] == [1, 1]:
            print(record)
```

## VCFWriter

To write a .vcf file you will need a `Schema` object, which is our representation of the vcf file header. This can either be passed into the VCFReaderContextManager constructor or a default minimal header will be written. A `Schema` is created on construction of the `VCFReader`.This is all done automatically if your header is well-formated.
 
Short example: open vcf file and write out all variant calls that were combined into a single line in the input vcf (multi alt calls)

```python
from vcfutils.parser import VCFReaderContextManager
from vcfutils.writer import VCFWriterContextManager

with VCFReaderContextManager('inFile.vcf') as vcf_reader:
    with VCFWriterContextManager('outFile.vcf', vcf_reader.header) as vcf_writer:
        for record in vcf_reader.read_records():
            if record.from_multi_alt:
                vcf_writer.write_record(record)
```

## VCF Adapters

VCF adapters are put BEFORE the read_record function and require a VCFReader object to work on. They can be put after each other in any order, but be aware that order matters.  In most cases you will want to sort the vcf at the end, i.e. in the outermost position.

The following adapters are available:
* `Sort`
* `SplitMNPs`
* `TrimPadding`
* `LeftAlign`
* `Unique` - after splitting of MNPs we can end up with two identical SNPs. This adapter removes the duplicate.
* `Standardise` - This is currently putting together (in the order of application) TrimPadding, SplitMNPs, LeftAlign, Sort, Unique.
* `Filter`

General usage:
Example: This will split MNPs when reading in a record:
```python
from echidna.vcfutils.parser import VCFReaderContextManager
from echidna.vcfutils import adapters

with VCFReaderContextManager('my_file.vcf') as vcf_reader:
    for record in adapters.SplitMNPs(vcf_reader).read_records():
        print(record)
```
Example: Multiple Adapters
```python
from echidna.vcfutils.parser import VCFReaderContextManager
from echidna.vcfutils import adapters

with VCFReaderContextManager('my_file.vcf') as vcf_reader:
    for record in adapters.Sort(adapters.SplitMNPs(vcf_reader)).read_records():
        print(record)
```
Example: LeftAlign requires a reference genome (fasta file)
```python
import pysam
from echidna.vcfutils.parser import VCFReaderContextManager
from echidna.vcfutils import adapters

with VCFReaderContextManager('my_file.vcf') as vcf_reader:
    for record in adapters.LeftAlign(vcf_reader, pysam.Fastafile('reference.fa')).read_records():
        print(record)
```
Example: Sort adapter sorts localy in a rolling window of 10 variants by default. To change the window width provide a second argument with the width to the Sort adapter.
```python
import pysam
from echidna.vcfutils.parser import VCFReaderContextManager
from echidna.vcfutils import adapters

with VCFReaderContextManager('my_file.vcf') as vcf_reader:
    for record in adapters.Sort(vcf_reader, window_size = 20).read_records():
        print(record)
```

Filter Adapter
--------------
A special type of adapter is a `Filter`.  The `Filter` adapter enables you to filter a vcf file for any characteristic, e.g. `snp_filter` allows you to extract all SNPs from a vcf file. A growing collection of filters can be found in 'vcfutils/filter.py'. Currently the collection includes e.g.
* `snp_filter`
* `ins_filter`
* `del_filter`
* `mnp_filter`
* `het_filter_for_sample(sample_name)`
* `hom_filter_for_sample(sample_name)`

Example use for a single filter (print all SNP's from a file):
```python
from echidna.vcfutils.parser import VCFReaderContextManager
from echidna.vcfutils import adapters
from echidna.vcfutils import filter

with VCFReaderContextManager(vcffile) as vcf_reader:
    for record in adapters.Filter(vcf_reader, filter.snp_filter).read_records():
        print(record)
```

Example use for a single filter - get all heterozygous variants for sample 'NA12878':
```python
from echidna.vcfutils.parser import VCFReaderContextManager
from echidna.vcfutils import adapters
from echidna.vcfutils import filter

with VCFReaderContextManager(vcffile) as vcf_reader:
    for record in adapters.Filter(vcf_reader, filter.het_filter_for_sample('NA12878')).read_records():
        print(record)
```

Filters can be combined using standard logical operators AND, OR, NOT as follows:

1)  Using operator NOT to extract all variants but SNP's
```python
from echidna.vcfutils.parser import VCFReaderContextManager
from echidna.vcfutils import adapters
from echidna.vcfutils import filter

with VCFReaderContextManager(vcffile) as vcf_reader:
    for record in adapters.Filter(vcf_reader, filter.filter_not(filter.snp_filter)).read_records():
        print(record)
```
2)  Use operator OR to extract all SNP's and MNP's
```python
from echidna.vcfutils.parser import VCFReaderContextManager
from echidna.vcfutils import adapters
from echidna.vcfutils import filter

with VCFReaderContextManager(vcffile) as vcf_reader:
    for record in adapters.Filter(
        vcf_reader,
        filter.filter_or(filter.snp_filter, filter.mnp_filter)
    ).read_records():
        print(record)
```
Another interesting example is filter for mendelian inconsistencies for trios. A ped file containing pedigree information about the family members has to be provided.
```python
with VCFReaderContextManager(vcffile) as vcf_reader:
    for record in adapters.Filter(
            vcf_reader,
            filter.filter_mendelian_inheritance_inconsistencies_from_trio('pedigree_file.ped')
    ).read_records():
        print(record)
```

If you want to filter on other properties, we strongly advise that you write a filter function into `vcfutils/filter.py`, so that it is reusable for other people. Each filter function needs to take a record as input and output a Boolean.

**Note:** `FilterWrapper` class is being deprecated.

## Comparison of vcfs
To compare two vcf files and obtain their intersection, right and left difference use the command line utility found and described in the [tools area](../../../tools).

The python interface of this utility is in the `comparison.py` module. Only variant information, i.e. chrom, pos, ref and alt are used for the comparison.  This is also the only information written to the output files.

The input files can be gzipped or raw vcf's. Optionally, an output filter string and output file dir can be passed into the constructor. Output filter string is then written into the FILTER field in the output files; PASS filter string is written otherwise. Output files `left_difference.vcf`, `right_difference.vcf` and `intersection.vcf` are written to the output dir provided or to the current working directory as default.  Headers and annotations for these files are inherited from the input files, the intersection inherits the 'left' file's header.

```python
from echidna.vcfutils.comparison import Comparison
comparator = Comparison(
    'left_hand_side_file.vcf',
    'right_hand_side_file.vcf',
    output_dir = optional_output_dir,
    filters = optional_filter_string
)
```
The output files have been saved at this point and you can access their location as `comparator.left_difference_filename` etc. You can also print out the stats about the comparison by printing the comparator: `print(comparator)`.

## TODO

- Support of the gVCF format.
- Serialisation into json.
- Ensure compatibility with global alliance standards.
- Split out variants into minimal units which can be glued together, e.g. with Haplotype ids. This is to be done in EchiDNA and here a support for this information has to be added.
