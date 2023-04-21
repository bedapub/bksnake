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
