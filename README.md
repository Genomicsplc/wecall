# GENOMICS plc - weCall
[!(Build status)(https://travis-ci.org/Genomicsplc/wecall.svg?branch=master)](https://travis-ci.org/Genomicsplc/wecall)

weCall is a **fast, accurate and simple** to use command line tool for **variant detection in Next Generation
Sequencing (NGS) data**.

It identifies genomic sites where genetic variants exist relative to the reference genome, and infers genotypes at these sites. 
It reads data from one or more BAM files simultaneously, each containing reads from one or more samples, and calls variants in these samples jointly. The output consists of a single VCF or gVCF file, and summarises the evidence for a variant site as a genotype call, a genotype quality, and genotype likelihoods for each sample.
weCall is fast to run and very easy to use in a single processing step. 

## Building

1. Ensure you have cmake, ncurses, boost and zlib installed.
2. Clone weCall repository. Withing the repository type ```make install```. This will add ```weCall``` and the unit test tools (```unittest```, ```iotest```)
to ```/usr/bin/local```. The installation directory can be set by modifying the ```PREFIX``` variable in the Makefile.
3. The user documentation can be created with ```make doc```. The user documentation will be generated in the folder ```build/doc/weCall-userguide.pdf``` of the repository.
4. weCall may be run directly from the command line. A quick summary of the parameters accepted by weCall can be obtained by invoking the help facility:
 ```weCall -h```.
For detailed information on how to use weCall, please refer to the user guide (```build/weCall-userguide.pdf```).

## Running weCall

weCall can be run directly from the command line. For an example of variant calling for chromosome 20: position 0 to 100000, run the following:
```bash
weCall --inputs sample.bam --refFile reference.fa --regions 20:0-100000
```
An `output.vcf` should appear in the working directory if it succeeded.

Multiple samples can be called simultanously by providing a comma-separeted list of bam files.
```bash
weCall --inputs sample1.bam,sample2.bam --refFile reference.fa --regions 20:0-100000
```

## User documentation

A quick summary of the parameters accepted by weCall may be obtained by invoking the help facility:
 ```weCall -h```.

For detailed information on how to use weCall, please refer to the user guide in ```build/weCall-userguide.pdf```. 
If it is not availble, it can be created by running ```make```.

## Testing Framework

weCall is tested using unit tests (C++) and acceptance test (python).

* To reformat the C++ code: ```./scripts/clang-reformat.sh```
* To run the C++ unit tests:  ```make test-unit```
* To run the full acceptance test suite you'll need python 3.x and a virtual environment. This will be automatically installed when you run the tests with ```make test-acceptance```.
