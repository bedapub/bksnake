# ----------------------------------------------------------------------------
"""
bksnake or vcsnake
here, we list all files necessary for the pipeline to have generated
prior to archive the data.
this also deletes all "empty" files, ie of size 0, in the "log" folder
"""
if config['pipeline'] == 'bksnake':
    rule processing_done:
        input:
            # Count data
            {TPM_COLLAPSED_GCT_REF},
            {TPM_COLLAPSED_GCT_ENS},
            #            
            # QC reports
            os.path.join(OD, 'multiqc_report.html'),
            os.path.join(OD, 'multiqc_report_ensembl.html'),
            os.path.join(OD_QC, 'refseq_log2tpm_pca.pdf'),
            os.path.join(OD_QC, 'bioQC.pdf'),
            #
            # Other optional output files (unmapped, bw, etc)
            expand(os.path.join(OD_CUTADAPT, '{sample}.report.txt'), sample=sample_ids),
            options.get_optional_output_files(sample_ids, config)
        output:
            os.path.join(OD_ANNO, 'processing.done')
        shell:
            """
            find {OD_LOG} -type f -empty -print -delete && touch {output[0]}
            """
elif config ['pipeline'] == 'vcsnake':
    rule processing_done:
        input:
            # QC reports
            os.path.join(OD, 'multiqc_report.html'),
            os.path.join(OD_VCF, 'allSamples_merged.vcf.gz'),
            os.path.join(OD_VCF, 'allSamples_merged.vcf.gz.tbi'),
            
            # Other optional output files (fastq)
            options.get_optional_output_files_vc(sample_ids, config)
        output:
            os.path.join(OD_ANNO, 'processing.done')
        shell:
            """
            find {OD_LOG} -type f -empty -print -delete && touch {output[0]}
            """
else:
    sys.exit('Variable \'pipeline\' not defined in input config yaml file.')
    


# ------------------------------------------------------------------
"""
Rule to remove output files that are not required
keep_fastq_files: True
keep_bam_files: True
"""

if config['keep_fastq_files'] != True:
    rule delete_fastq:
        input:
            os.path.join(OD_FASTQ, '{name}.fastq.gz'),
            os.path.join(OD, 'multiqc_report.html')
        output:
            temp(os.path.join(OD_FASTQ, '{name}.fastq.gz_to_delete'))
        threads: 1
        resources:
            mem_mb=1000
        shell:
            """
            rm -f {input[0]} && touch {output}
            """


if config['keep_bam_files'] != True:
    rule delete_bam:
        input:
            os.path.join(OD_BAM, '{name}.bam'),
            os.path.join(OD_BAM, '{name}.bam.bai'),
            os.path.join(OD_BW, '{name}.bw.done'),
            os.path.join(OD, 'multiqc_report.html')
        output:
            temp(os.path.join(OD_BAM, '{name}.bam_to_delete')),
            temp(os.path.join(OD_BAM, '{name}.bam.bai_to_delete'))
        threads: 1
        resources:
            mem_mb=1000
        shell:
            """
            rm -f {input[0]} && touch {output[0]}
            rm -f {input[1]} && touch {output[1]}
            """
