# bksnake
Public version of bksnake - biokit snakemake - bulk RNASeq Snakemake workflow



## Introduction

Snakemake ([Moelder et al., 2021](https://f1000research.com/articles/10-33/v1)) implements a bulk RNASeq data analysis workflow using STAR aligner ([Dobin et al., 2012](https://academic.oup.com/bioinformatics/article/29/1/15/272537)) for read mapping and FeatureCounts from the Subread package ([Liao et al., 2014](https://pubmed.ncbi.nlm.nih.gov/24227677/)) for gene quantification. Reference genomes with RefSeq and Ensembl gene annotations are available for several species such as hg38, chm13, mm10, mm39, rn6, rn7, mfa5, mfa6, ss11, and oc2. The generation of these reference genomes and annotation files is documented in a separate repository that is currently under construction. Data quality and RNASeq metrics are determined using FastQC [Andrews et al.](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)), MultiQC ([Ewels et al., 2015](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507)), and Picard tools ([Broad Institute](http://broadinstitute.github.io/picard/)). In addition, diagnostic plots for data quality assessment such as BioQC tissue heterogeneity ([Zhang et al., 2017](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3661-2)) or Principal Component Analysis are provided in an HTML report that is also currently under construction. Optionally, genome coverage files (BigWig) and read alignment files (BAM/CRAM) can be generated as well. Input read trimming with Cutadapt ([Martin, 2010](https://cutadapt.readthedocs.io/en/stable)) and generation of unmapped reads are also available. The pipeline can be launched via a helper tool, run.py, or directly with Snakemake for users familiar with the workflow tool. All parameters for the pipeline are specified within a configuration yaml file or explicitly on the command line when using run.py. All input data, i.e. input fastq files, a human-readable tab-delimited file describing the samples, as well as the reference genome and STAR index files, must be available to the pipeline in a local data folder. To run the pipeline, Snakemake and [Singularity](https://sylabs.io/docs/) must be installed and pre-configured. All software tools used by the pipeline are pulled from public Singularity or Docker image repositories. It is recommended to run the pipeline on a high-performance cluster environment.


### Overview of the analysis workflow

Directed acyclic graph of a representative analysis workflow.

<div
  <figure>
    <img src="./resources/dag.png" alt="DAG of the workflow." style="width:50%">
    <figcaption>DAG of the workflow.</figcaption>
  </figure>
</div>


## Requirements

This workflow requires the following tools in the system path

1. [Singularity](https://docs.sylabs.io/guides/main/user-guide/)
2. [Snakemake](https://snakemake.readthedocs.io/en/stable/)

In order to pull Singularity images from GitHub package registry one needs to specify username and GitHub read package token:

```bash

export SINGULARITY_DOCKER_USERNAME=<username>
export SINGULARITY_DOCKER_PASSWORD=<github read package token>

```

## Preparation

### Reference genome (TBD)

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

It is possible to add more columns, for example, to describe additional experiment parameters or specimen information. But they are not used by the pipeline further.

#### Example

| #ID        | GROUP       | FASTQ1                     | FASTQ2                     | Raw                       | Organism |
|------------|-------------|----------------------------|----------------------------|---------------------------|----------|
| GSM5362225 | HOXB9-KO    | 10k_SRR14749833_1.fastq.gz | 10k_SRR14749833_2.fastq.gz | resources/test-data/fastq | Human    |
| GSM5362224 | HOXB9-T133A | 10k_SRR14749832_1.fastq.gz | 10k_SRR14749832_2.fastq.gz | resources/test-data/fastq | Human    |
| GSM5362223 | HOXB9       | 10k_SRR14749831_1.fastq.gz | 10k_SRR14749831_2.fastq.gz | resources/test-data/fastq | Human    |

### Configuration

The workflow requires several parameters to be configured, most of which can be set through a yaml configuration file. 
A template file named `config.yaml` is provided in the config directory. 
Note that many of these parameters can also be specified through the wrapper script `run.py`, as explained in the next section. 
Parameters specified on the command line through the wrapper script will overwrite parameters set in the configuration file.

To learn about all possible parameters, execute:

```bash

python run.py --help

```

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

Here is the folder structure of a typical workflow run

```

├── annot                            genome annotation files
├── bam                              read mappings in BAM format (optional)
├── bw                               read coverage BigWig files (optional)
├── config.yaml                      workflow configuration file
├── cram                             read mappings in CRAM format (optional)
├── cutadapt                         Cudapapt output (optional)
├── fastq                            copy of input reads (optional)
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
├── samples.txt                      sample metadata in tab-delimited file
└── unmapped                         unmapped reads in FASTQ file format (optional)

```

Explanations

### "annot" folder

Contains copies of reference genome and annotations files (i.e. `FASTA`, `GTF`, etc)
In addition, sample metadata in tab-delimited text file, `phenoData.meta` and in [`cls`](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#CLS:_Continuous_.28e.g_time-series_or_gene_profile.29_file_format_.28.2A.cls.29) format `phenoData.cls`.

### "bam" folder (optional)

Contains all aligned reads in [`BAM`](https://en.wikipedia.org/wiki/BAM_(file_format)) file. Only generated if the parameter `keep_bam_files` is `True` in the pipeline configuration.

### "bw" folder (optional)

Read coverage files in [BigWig](http://genome.ucsc.edu/goldenPath/help/bigWig.html) format. May be used for graphical visualisation of the genome coverage by external tools such as [JBrowse](https://jbrowse.org/jb2/docs/user_guide/) or [IGV](https://software.broadinstitute.org/software/igv/userguide). Only generated if the parameter `generate_bw_files` is `True` in the pipeline configuration.

### "cutadapt" folder (optional)

Output files from the sequencing read trimming tool [`Cutadapt`](https://cutadapt.readthedocs.io/en/stable/) ([Martin 2010](https://journal.embnet.org/index.php/embnetjournal/article/view/200)). Only generated if the parameter `cutadapt: run` is `True` in the pipeline configuration.

### "config.yaml" file

All configuration parameters stored in a single [`yaml`](https://en.wikipedia.org/wiki/YAML) file.

### "cram" folder (optional)

Contains all aligned reads in [`CRAM`](https://en.wikipedia.org/wiki/CRAM_(file_format)) file. Only generated if the parameter `generate_cram_files` is `True` in the pipeline configuration.

### "fastq" folder (optional)

Contains a copy of the input [`FASTQ`](https://en.wikipedia.org/wiki/FASTQ_format) files. Only generated if the parameter `keep_fastq_files` is `True` in the pipeline configuration.

### "fastqc" folder

Contains all output files from the read quality control tool [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

### "fc" folder

Intermediate output files, summary and compressed counts file, from the gene quantification step by using [`FeatureCounts`](https://subread.sourceforge.net/featureCounts.html) tool from the `Subreads` package. Parameters for `FeatureCounts` can be specified via the input configuration file.

### "gct" folder

Contains _RefSeq_ annotated gene counts, normalized counts, _log2_-transformed counts in [`GCT`](https://software.broadinstitute.org/software/igv/GCT) file format.

- `refseq_count.gct`: RefSeq gene counts
- `refseq_count_collapsed.gct`: RefSeq gene counts collapsed to human orthologous gene symbols (using `resources/geneids.chip`)
- `refseq_tpm.gct`: normalized RefSeq transcript per million mapped reads (tpm)
- `refseq_tpm_collapsed.gct`: human orthologs of normalized counts (tpm)
- `refseq_log2tpm.gct`: log2-transformed normalized RefSeq transcript per million mapped reads (tpm)

### "gct-ens" folder

Contains _Ensembl_ annotated gene counts, normalized counts, _log2_-transformed counts in [`GCT`](https://software.broadinstitute.org/software/igv/GCT) file format.

- `ensembl_count.gct`: Ensembl gene counts
- `ensembl_count_collapsed.gct`: Ensembl counts collapsed to human orthologous gene symbols (using `resources/ENSEMBLGENES.chip`)
- `ensembl_tpm.gct`: normalized Ensembl transcript per million mapped reads (tpm)
- `ensembl.gct`: human orthologs of normalized counts (tpm)
- `ensembl_log2tpm.gct`: log2-transformed normalized RefSeq transcript per million mapped reads (tpm)

### "log" folder

Contains several log files from analysis tools used by the pipeline. Mainly used for debugging purposes.

### "multiqc_data" folder

Output files from the [`MultiQC`](https://multiqc.info/docs/) tool with _RefSeq_ gene annotations (e.g. for Picard RNASeq metrics).

### "multiqc_data_ensembl" folder

Output files from the [`MultiQC`](https://multiqc.info/docs/) tool with _Ensembl_ gene annotations (e.g. for Picard RNASeq metrics).

### "multiqc_report_ensembl.html" file

HTML summary report from the [`MultiQC`](https://multiqc.info/docs/) tool based on _Ensembl_ genome annotations.

### "multiqc_report.html" file

HTML summary report from the [`MultiQC`](https://multiqc.info/docs/) tool based on _RefSeq_ genome annotations.

### "qc" folder

Plots for quality control purposes.

- `bioQC.pdf`: Heatmap representationg of the _BioQC_ enrichment scores for detecting detecting such tissue heterogeneity ([Zhang et al. 2017](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3661-2))
- `bioQC_thr2.txt`: _BioQC_ enrichment scores above threshold 2
- `bioQC.txt`: All _BioQC_ enrichment scores 
- `refseq_log2tpm_pca.pdf`: Plot of the main components from the principal component analysis on the basis of the _log2_-transformed normalized _RefSeq_ gene counts (log2tpm).
- `refseq_log2tpm_pca.txt`: Coordinates of the components from the principal component analysis on the basis of the _log2_-transformed normalized _RefSeq_ gene counts (log2tpm).

### "rulegraph.pdf" and "rulegraph.png" files

Graphical representation of the entire workflow. A directed acyclic graph, DAG, generated by [`Snakemake`](https://snakemake.readthedocs.io/en/stable/)

### "samples.txt" file

Tab-delimited, human readable text file containing study and sample metadata in tabular form (one sample per line).

### "unmapped" folder (optional)

Contains unmapped sequencing reads in [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files. Only generated if the parameter `generate_unmapped` is `True` in the pipeline configuration.
