"""
functions to handle options
e.g. optional output files
"""
import os
import sys
import pandas as pd
from snakemake.io import expand

# ------------------------------------------------------------------------------
# Declare all required output files in a list object. This will be given to the rule 'all'
# For bksnake / bulk rnaseq pipeline
def get_optional_output_files(sample_ids, config):

    OD = config['outdir']
    OD_FASTQ = os.path.join(OD, 'fastq')
    OD_BAM = os.path.join(OD, 'bam')
    OD_CRAM = os.path.join(OD, 'cram')
    OD_BW = os.path.join(OD, 'bw')

    optional_output_files = []
    if config['keep_fastq_files'] == True:
        if config['library']['type'] == 'paired-end':
            optional_output_files.append(expand(os.path.join(OD_FASTQ, '{sample}_1.fastq.gz'), sample=sample_ids))
            optional_output_files.append(expand(os.path.join(OD_FASTQ, '{sample}_2.fastq.gz'), sample=sample_ids))
        else:
            optional_output_files.append(expand(os.path.join(OD_FASTQ, '{sample}.fastq.gz'), sample=sample_ids))
    else:
        if config['library']['type'] == 'paired-end':
            optional_output_files.append(expand(os.path.join(OD_FASTQ, '{sample}_1.fastq.gz_to_delete'), sample=sample_ids))
            optional_output_files.append(expand(os.path.join(OD_FASTQ, '{sample}_2.fastq.gz_to_delete'), sample=sample_ids))
        else:
            optional_output_files.append(expand(os.path.join(OD_FASTQ, '{sample}.fastq.gz_to_delete'), sample=sample_ids))


    if config['keep_bam_files'] == True:
        optional_output_files.append(expand(os.path.join(OD_BAM, '{sample}.bam'), sample=sample_ids))
        optional_output_files.append(expand(os.path.join(OD_BAM, '{sample}.bam.bai'), sample=sample_ids))
    else:
        optional_output_files.append(expand(os.path.join(OD_BAM, '{sample}.bam_to_delete'), sample=sample_ids))
        optional_output_files.append(expand(os.path.join(OD_BAM, '{sample}.bam.bai_to_delete'), sample=sample_ids))


    if config['generate_cram_files'] == True:
        optional_output_files.append(expand(os.path.join(OD_CRAM, '{sample}.cram'), sample=sample_ids))
        optional_output_files.append(expand(os.path.join(OD_CRAM, '{sample}.cram.crai'), sample=sample_ids))


    if config['generate_bw_files'] == True:
        optional_output_files.append(expand(os.path.join(OD_BW, '{sample}.bw'), sample=sample_ids))


    if config['generate_unmapped'] == True:
        if config['library']['type'] == 'paired-end':
            optional_output_files.append(expand(os.path.join(OD, 'unmapped', '{sample}_1.fastq.gz'), sample=sample_ids))
            optional_output_files.append(expand(os.path.join(OD, 'unmapped', '{sample}_2.fastq.gz'), sample=sample_ids))
        else:
            optional_output_files.append(expand(os.path.join(OD, 'unmapped', '{sample}_1.fastq.gz'), sample=sample_ids))

    return optional_output_files
#print('List of all optional output files')
#from pprint import pprint
#pprint(optional_output_files)



# ------------------------------------------------------------------------------
# Declare all required output files in a list object. This will be given to the rule 'all'
# For vcsnake / variant calling pipeline
def get_optional_output_files_vc(sample_ids, config):

    OD = config['outdir']
    OD_FASTQ = os.path.join(OD, 'fastq')
    OD_BAM = os.path.join(OD, 'bam')
    OD_CRAM = os.path.join(OD, 'cram')
    #OD_BW = os.path.join(OD, 'bw')

    optional_output_files = []
    if config['keep_fastq_files'] == True:
        if config['library']['type'] == 'paired-end':
            optional_output_files.append(expand(os.path.join(OD_FASTQ, '{sample}_1.fastq.gz'), sample=sample_ids))
            optional_output_files.append(expand(os.path.join(OD_FASTQ, '{sample}_2.fastq.gz'), sample=sample_ids))
        else:
            optional_output_files.append(expand(os.path.join(OD_FASTQ, '{sample}.fastq.gz'), sample=sample_ids))
    else:
        if config['library']['type'] == 'paired-end':
            optional_output_files.append(expand(os.path.join(OD_FASTQ, '{sample}_1.fastq.gz_to_delete'), sample=sample_ids))
            optional_output_files.append(expand(os.path.join(OD_FASTQ, '{sample}_2.fastq.gz_to_delete'), sample=sample_ids))
        else:
            optional_output_files.append(expand(os.path.join(OD_FASTQ, '{sample}.fastq.gz_to_delete'), sample=sample_ids))


    if config['keep_bam_files'] == True:
        optional_output_files.append(expand(os.path.join(OD_BAM, '{sample}.bam'), sample=sample_ids))
        optional_output_files.append(expand(os.path.join(OD_BAM, '{sample}.bam.bai'), sample=sample_ids))
    else:
        optional_output_files.append(expand(os.path.join(OD_BAM, '{sample}.bam_to_delete'), sample=sample_ids))
        optional_output_files.append(expand(os.path.join(OD_BAM, '{sample}.bam.bai_to_delete'), sample=sample_ids))


    if config['generate_cram_files'] == True:
        optional_output_files.append(expand(os.path.join(OD_CRAM, '{sample}.cram'), sample=sample_ids))
        optional_output_files.append(expand(os.path.join(OD_CRAM, '{sample}.cram.crai'), sample=sample_ids))


    #if config['generate_bw_files'] == True:
    #    optional_output_files.append(expand(os.path.join(OD_BW, '{sample}.bw'), sample=sample_ids))


    #if config['generate_unmapped'] == True:
    #    if config['library']['type'] == 'paired-end':
    #        optional_output_files.append(expand(os.path.join(OD, 'unmapped', '{sample}_1.fastq.gz'), sample=sample_ids))
    #        optional_output_files.append(expand(os.path.join(OD, 'unmapped', '{sample}_2.fastq.gz'), sample=sample_ids))
    #    else:
    #        optional_output_files.append(expand(os.path.join(OD, 'unmapped', '{sample}_1.fastq.gz'), sample=sample_ids))

    return optional_output_files
#print('List of all optional output files')
#from pprint import pprint
#pprint(optional_output_files)
