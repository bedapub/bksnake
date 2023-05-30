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

_**TBD** This section is still under construction..._

Structure of the genomes directory
```
/path/to/genome/root/directory/hg38/
├── fasta
├── gff3
│   ├── ensembl
│   └── refseq
├── gtf
│   ├── ensembl
│   └── refseq
├── star_2.7.10b
```

In each "sub" directory, named with genome version (e.g. hg38), there are subfolders for `fasta`, `gff3`, `gtf3` and `STAR` index files.


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
Most of these parameters can be configured by a "yaml" configuration file.
A template file is given by `config/config.yaml`.
Note that many of these parameters can also be specified via the wrapper script `run.py`, see next section.
Parameters specified on the command line via the wrapper script will overwrite parameters in the configuration file.

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

python run.py \
    --snakemake-path="ml purge && ml snakemake && snakemake" \
    --genome-dir /data/genomes \
    --outdir output \
    --species hg38 \
    --jobs 100 

```

## Output

Here is the folder structure of a typical workflow run

```bash

├── annot                            genome annotation files
├── config.yaml                      workflow configuration file
├── cram                             read mappings in CRAM format
├── fastqc                           FASTQC output files
├── fc                               FeatureCounts output files
├── gct                              gene counts and normalized gene counts in GCT file format for RefSeq annotations
├── gct-ens                          gene counts and normalized gene counts in GCT file format for Ensembl annotations
├── log                              log and output files from the tools used
├── multiqc_data                     MultiQC data files for RefSeq annotations
├── multiqc_data_ensembl             MultiQC data files for Ensembl annotations
├── multiqc_report_ensembl.html      MultiQC report for Ensembl annotations
├── multiqc_report.html              MultiQC report for RefSeq annotations
├── qc                               some QC plots, e.g. PCA or BioQC
├── rulegraph.pdf                    workflow DAG in pdf format
├── rulegraph.png                    workflow DAG in png format
└── samples.txt                      sample metadata in tab-delimited file

```

Explanations of each file and folder.

### "annot" folder

Contains copies of reference genome and annotations files.
In addition, sample metadata in tab-delimited text file, `phenoData.meta` and in [`cls`](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#CLS:_Continuous_.28e.g_time-series_or_gene_profile.29_file_format_.28.2A.cls.29) format `phenoData.cls`.

### "config.yaml" file

All configuration parameters on one `yaml` file.

### "cram" folder

Contains all aligned reads in [`CRAM`](https://en.wikipedia.org/wiki/CRAM_(file_format)) file.

### "fastqc" folder

Contains all output files from the read quality control tool [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

### "qc" folder

Contains some plots for quality control purposes.

- `bioQC.pdf`: Heatmap representationg of the BioQC enrichment scores for detecting detecting such tissue heterogeneity ([Zhang et al. 2017](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3661-2))
- `bioQC_thr2.txt`: BioQC enrichment scores above threshold 2
- `bioQC.txt`: All BioQC enrichment scores 
- `refseq_log2tpm_pca.pdf`: Plot of the main components from the principal component analysis on the basis of the log2-transformed normalized RefSeq gene counts (log2tpm).
- `refseq_log2tpm_pca.txt`: Coordinates of the components from the principal component analysis on the basis of the log2-transformed normalized RefSeq gene counts (log2tpm).


### "gct" folder

Contains RefSeq annotated gene counts, normalized counts, log2-transformed counts in [`GCT`](https://software.broadinstitute.org/software/igv/GCT) file format.

- `refseq_count.gct`: RefSeq gene counts
- `refseq_count_collapsed.gct`: RefSeq gene counts collapsed to human orthologous gene symbols (using `resources/geneids.chip`)
- `refseq_tpm.gct`: normalized RefSeq transcript per million mapped reads (tpm)
- `refseq_tpm_collapsed.gct`: human orthologs of normalized counts (tpm)
- `refseq_log2tpm.gct`: log2-transformed normalized RefSeq transcript per million mapped reads (tpm)

### "gct-ens" folder

Contains Ensembl annotated gene counts, normalized counts, log2-transformed counts in [`GCT`](https://software.broadinstitute.org/software/igv/GCT) file format.

- `ensembl_count.gct`: Ensembl gene counts
- `ensembl_count_collapsed.gct`: Ensembl counts collapsed to human orthologous gene symbols (using `resources/ENSEMBLGENES.chip`)
- `ensembl_tpm.gct`: normalized Ensembl transcript per million mapped reads (tpm)
- `ensembl.gct`: human orthologs of normalized counts (tpm)
- `ensembl_log2tpm.gct`: log2-transformed normalized RefSeq transcript per million mapped reads (tpm)

### "log" folder

Contains several log file from various tools used by the pipeline. Mainly used for debugging purposes.

### "multiqc_data" folder

Contains

### "multiqc_data_ensembl" folder
### "multiqc_report_ensembl.html" file
### "multiqc_report.html" file
### "qc" folder
### "rulegraph.pdf" file
### "rulegraph.png" file
### "samples.txt" file

