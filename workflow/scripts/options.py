"""
functions to handle options
e.g. optional output files
"""
import os
import sys
import pandas as pd
from snakemake.io import expand


# ----------------------------------------------------------------------------
"""
bksnake or vcsnake. Currently, only "bksnake" is supported.
here, we list all files necessary for the pipeline to have generated
prior to archive the data.
this also deletes all "empty" files, ie of size 0, in the "log" folder
"""
def get_all_output_files(config, sample_ids, dbs):

    OD = config['results']
    OD_QC = os.path.join(OD, 'qc')
    OD_GCT = os.path.join(OD, 'gct')
    OD_VCF = os.path.join(OD, 'vcf')

    all_output_files = []

    if config['pipeline'] == 'bksnake':
        # Count data
        all_output_files = expand(os.path.join(OD_GCT, '{db}_tpm_collapsed.gct'), db=dbs)
        
        # QC reports
        all_output_files += expand(os.path.join(OD, '{db}_multiqc_report.html'), db=dbs)
        all_output_files += expand(os.path.join(OD_QC, '{db}_log2tpm_pca.png'), db=dbs)
        all_output_files += expand(os.path.join(OD_QC, '{db}_bioQC.png'), db=dbs)
           
        # Other optional output files (unmapped, bw, etc)
        #all_output_files += expand(os.path.join(OD_CUTADAPT, '{sample}.report.txt'), sample=sample_ids)
        all_output_files += get_optional_output_files(sample_ids, config, dbs)
        
    elif config ['pipeline'] == 'vcsnake':
        # QC reports
        all_output_files.append(os.path.join(OD, 'multiqc_report.html'))
        all_output_files.append(os.path.join(OD_VCF, 'allSamples_merged.vcf.gz'))
        all_output_files.append(os.path.join(OD_VCF, 'allSamples_merged.vcf.gz.tbi'))

        # Other optional output files (fastq)
        all_output_files += get_optional_output_files_vc(sample_ids, config, dbs)
        
    else:
        sys.exit('Variable \'pipeline\' not defined in input config yaml file.')

    return all_output_files
        

# ------------------------------------------------------------------------------
# Declare all required output files in a list object. This will be given to the rule 'all'
# For bksnake / bulk rnaseq pipeline
def get_optional_output_files(sample_ids, config, dbs):

    OD = config['results']
    OD_BW = os.path.join(OD, 'bw')
    OD_BAM = os.path.join(OD, 'bam')
    OD_CRAM = os.path.join(OD, 'cram')
    OD_FASTQ = os.path.join(OD, 'fastq')
    OD_METRICS = os.path.join(OD, 'metrics')

    optional_output_files = []
    
    if config['keep_fastq_files'] == True:
        if config['library']['type'] == 'paired-end':
            optional_output_files += expand(os.path.join(OD_FASTQ, '{sample}_1.fastq.gz'), sample=sample_ids)
            optional_output_files += expand(os.path.join(OD_FASTQ, '{sample}_2.fastq.gz'), sample=sample_ids)
        else:
            optional_output_files += expand(os.path.join(OD_FASTQ, '{sample}.fastq.gz'), sample=sample_ids)

    if config['keep_bam_files'] == True:
        optional_output_files += expand(os.path.join(OD_BAM, '{sample}.bam'), sample=sample_ids)
        optional_output_files += expand(os.path.join(OD_BAM, '{sample}.bam.bai'), sample=sample_ids)

    if config['generate_cram_files'] == True:
        optional_output_files += expand(os.path.join(OD_CRAM, '{sample}.cram'), sample=sample_ids)
        optional_output_files += expand(os.path.join(OD_CRAM, '{sample}.cram.crai'), sample=sample_ids)

    if config['generate_bw_files'] == True:
        optional_output_files += expand(os.path.join(OD_BW, '{sample}.bw'), sample=sample_ids)

    if config['generate_unmapped'] == True:
        if config['library']['type'] == 'paired-end':
            optional_output_files += expand(os.path.join(OD, 'unmapped', '{sample}_1.fastq.gz'), sample=sample_ids)
            optional_output_files += expand(os.path.join(OD, 'unmapped', '{sample}_2.fastq.gz'), sample=sample_ids)
        else:
            optional_output_files += expand(os.path.join(OD, 'unmapped', '{sample}_1.fastq.gz'), sample=sample_ids)

    if config['junction_annotation'] == True:
        optional_output_files += expand(os.path.join(OD_METRICS, '{sample}.{db}.junction_annotation.log'), sample=sample_ids, db=dbs)
        optional_output_files += expand(os.path.join(OD_METRICS, '{sample}.{db}.splice_events.pdf'), sample=sample_ids, db=dbs)
        optional_output_files += expand(os.path.join(OD_METRICS, '{sample}.{db}.splice_junction.pdf'), sample=sample_ids, db=dbs)
        optional_output_files += expand(os.path.join(OD_METRICS, '{sample}.{db}.junction_plot.r'), sample=sample_ids, db=dbs)
        optional_output_files += expand(os.path.join(OD_METRICS, '{sample}.{db}.junction.xls'), sample=sample_ids, db=dbs)
        optional_output_files += expand(os.path.join(OD_METRICS, '{sample}.{db}.junction.bed'), sample=sample_ids, db=dbs)
        optional_output_files += expand(os.path.join(OD_METRICS, '{sample}.{db}.junction.Interact.bed'), sample=sample_ids, db=dbs)
        optional_output_files += expand(os.path.join(OD_METRICS, '{sample}.{db}.junctionSaturation_plot.pdf'), sample=sample_ids, db=dbs)
        optional_output_files += expand(os.path.join(OD_METRICS, '{sample}.{db}.junctionSaturation_plot.r'), sample=sample_ids, db=dbs)

# Currently, fingerprinting works only for the hg38/GRCh38p14 genome
    # In order to use the T2T_CHM13v2_0 genome, a new haplotype.map file needs to be created.
    if config['organism'] == 'Homo sapiens' and (config['species'] == 'GRCh38p14' or config['species'] == 'hg38') and 'crosscheck_fingerprints' in config:
        if config['crosscheck_fingerprints'] == True:
            optional_output_files += [os.path.join(OD_METRICS, 'crosscheck_metrics')]

    return optional_output_files


# ------------------------------------------------------------------------------
# Declare all required output files in a list object. This will be given to the rule 'all'
# For vcsnake / variant calling pipeline
def get_optional_output_files_vc(sample_ids, config):

    OD = config['results']
    OD_BAM = os.path.join(OD, 'bam')
    OD_CRAM = os.path.join(OD, 'cram')
    OD_FASTQ = os.path.join(OD, 'fastq')

    optional_output_files = []
    
    if config['keep_fastq_files'] == True:
        if config['library']['type'] == 'paired-end':
            optional_output_files += expand(os.path.join(OD_FASTQ, '{sample}_1.fastq.gz'), sample=sample_ids)
            optional_output_files += expand(os.path.join(OD_FASTQ, '{sample}_2.fastq.gz'), sample=sample_ids)
        else:
            optional_output_files += expand(os.path.join(OD_FASTQ, '{sample}.fastq.gz'), sample=sample_ids)

    if config['keep_bam_files'] == True:
        optional_output_files += expand(os.path.join(OD_BAM, '{sample}.bam'), sample=sample_ids)
        optional_output_files += expand(os.path.join(OD_BAM, '{sample}.bam.bai'), sample=sample_ids)

    if config['generate_cram_files'] == True:
        optional_output_files += expand(os.path.join(OD_CRAM, '{sample}.cram'), sample=sample_ids)
        optional_output_files += expand(os.path.join(OD_CRAM, '{sample}.cram.crai'), sample=sample_ids)

    return optional_output_files
