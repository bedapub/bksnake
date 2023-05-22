# bksnake
Public version of bksnake - biokit snakemake - bulk RNASeq Snakemake workflow



## Introduction

### Overview of the analysis workflow

<div
  <figure>
    <img src="./resources/dag.png" alt="DAG of the workflow." style="width:50%">
    <figcaption>DAG of the workflow.</figcaption>
  </figure>
</div>


## Installation

In order to pull Singularity images from GitHub package registry one needs to specify username and GitHub read package token:

```bash

export SINGULARITY_DOCKER_USERNAME=<username>
export SINGULARITY_DOCKER_PASSWORD=<github read package token>

```

## Preparation

### Tools

This workflow requires the following two tools in the system path

- singularity
- snakemake
- optional: cluster queue configuration

### Reference genome

Download reference genome annotation files.

### Metadata

- create metadata file with sample and fastq file information
- required columns are
- #ID: unique sample label
- GROUP: name of sample condition
- FASTQ1: name of forward read fastq file (mate 1, R1)
- FASTQ2: name of reverse read fastq file (mate 2, R2)
- Raw: path to the folder containing the input fastq files
- Organism: Human, Rat, Mouse, or Pig
- use test data

#### Example

| #ID        | GROUP       | FASTQ1                     | FASTQ2                     | Raw                       | Organism |
|------------|-------------|----------------------------|----------------------------|---------------------------|----------|
| GSM5362225 | HOXB9-KO    | 10k_SRR14749833_1.fastq.gz | 10k_SRR14749833_2.fastq.gz | resources/test-data/fastq | Human    |
| GSM5362224 | HOXB9-T133A | 10k_SRR14749832_1.fastq.gz | 10k_SRR14749832_2.fastq.gz | resources/test-data/fastq | Human    |
| GSM5362223 | HOXB9       | 10k_SRR14749831_1.fastq.gz | 10k_SRR14749831_2.fastq.gz | resources/test-data/fastq | Human    |

### Configuration

Edit configuration file `config/config.yaml`

- genome directory
- singularity directory
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
 
