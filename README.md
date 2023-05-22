# bksnake
Public version of bksnake - biokit snakemake - bulk RNASeq Snakemake workflow



## Introduction

### Overview of the analysis workflow

Directed acyclic graph of a representative analysis workflow.

<div
  <figure>
    <img src="./resources/dag.png" alt="DAG of the workflow." style="width:50%">
    <figcaption>DAG of the workflow.</figcaption>
  </figure>
</div>


## Installation

This workflow requires the following tools in the system path

1. [Singularity](https://docs.sylabs.io/guides/main/user-guide/)
2. [Snakemake](https://snakemake.readthedocs.io/en/stable/)

In order to pull Singularity images from GitHub package registry one needs to specify username and GitHub read package token:

```bash

export SINGULARITY_DOCKER_USERNAME=<username>
export SINGULARITY_DOCKER_PASSWORD=<github read package token>

```

## Preparation

### Reference genome

Download reference genome annotation files into a genome "root" and "sub" directories.
Specify these directories also in the pipeline configuration file (see below).
TBD

### Metadata

The main input to the workflow is a tab-delimited text file containing metadata information about all samples. 
There is a header line and one sample per line, organized by several columns as follows:

- `#ID`: unique sample label
- `GROUP`: name of sample condition
- `FASTQ1`: name of forward read fastq file (mate 1, R1)
- `FASTQ2`: name of reverse read fastq file (mate 2, R2)
- `Raw`: path to the folder containing the input fastq files
- `Organism`: Human, Rat, Mouse, or Pig

Note that columns `#ID` and `GROUP` may not contain white-spaces and other special characters.

#### Example

| #ID        | GROUP       | FASTQ1                     | FASTQ2                     | Raw                       | Organism |
|------------|-------------|----------------------------|----------------------------|---------------------------|----------|
| GSM5362225 | HOXB9-KO    | 10k_SRR14749833_1.fastq.gz | 10k_SRR14749833_2.fastq.gz | resources/test-data/fastq | Human    |
| GSM5362224 | HOXB9-T133A | 10k_SRR14749832_1.fastq.gz | 10k_SRR14749832_2.fastq.gz | resources/test-data/fastq | Human    |
| GSM5362223 | HOXB9       | 10k_SRR14749831_1.fastq.gz | 10k_SRR14749831_2.fastq.gz | resources/test-data/fastq | Human    |

### Configuration

The workflow requires several parameters to be configured.
Most of these parameters can be configured by a yaml configuration file.
A template file is given by `config/config.yaml`.
Note that many of these parameters can also be specified via the wrapper script `run.py`, see next section.
Parameters specified on the command line via the wrapper script over write parameters in the configuration file.

It is important to specify the following parameters
- genome root and sub- directory
- singularity images directory
- snakemake path
- sequencing library (single-end or paired-end, stranded or unstranded)

## Usage

Run on a cluster with LSF scheduler, up to 100 jobs in parallel

```bash

python run.py --outdir output --species hg38 --jobs 100

```

Run locally, using up to 12 cores


```bash

python run.py --outdir output --species hg38 --cores 12

```
 
Here, an example where several pipeline parameters are specified directly via the command line
```bash

python run.py --snakemake-path="ml purge && ml snakemake && snakemake" --genome-dir /data/genomes --outdir output --species hg38 --jobs 100 

```