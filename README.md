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

### Configuration

Edit configuration file `config/config.yaml`

- genome directory
- singularity directory
- snakemake path
- sequencing library

### Metadata

- create metadata file with sample and fastq file information
- use test data

## Usage

Run on a cluster with LSF scheduler, up to 100 jobs in parallel

```bash

python run.py --outdir output --species hg38 --jobs 100

```

Run locally, using up to 12 cores


```bash

python run.py --outdir output --species hg38 --cores 12

```
 
