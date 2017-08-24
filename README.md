rna-tpm-calculator
==================

A TPM calculating workflow for multiple samples in RNA-seq

This workflow depends:

- Bowtie2
- FastQC
- pigz
- PRINSEQ
- R (with docopt and tidyverse)
- RSEM
- SAMtools

Usage
-----

```sh
$ ./workflow.sh --help
Usage:  workflow.sh [ -h | --help | -v | --version ]
        workflow.sh run [--thread <int>] [--no-qc] [--ref <fasta>]
                        [--ref-gene-map <txt>] [--input <dir>]
                        [--output <dir>]

Description:
  TPM calculating workflow for RNA-seq

Arguments:
  -h, --help              Print usage
  -v, --version           Print version information
  --thread <int>          Limit multithreading
  --no-qc                 Skip QC checks with FastQC
  --ref <fasta>           Pass a gzip-compressed reference FASTA file
                          [default: ./ref/GRCm38_p5.fasta.gz]
  --ref-gene-map <txt>    Pass a gzip-compressed transcript-to-gene-map file
                          [default: ./ref/GRCm38_p5_gene_map.txt.gz]
  --input <dir>           Pass an path to input directory
                          [default: ./input]
  --output <dir>          Pass an path to output directoryr
                          [default: ./output]
  run                     Run the workflow
```

Docker image
------------

```sh
# build using docker
$ docker image build -t tpm ./build

# build using docker-compose
$ cp example-docker-compose.yml docker-compose.yml
$ docker-compose build
```
