# ------------------------------------------------------------------------------
# Picard mRNA metrics only for paired-end reads
rule picard_ensembl:
    input:
        os.path.join(OD_BAM,'{sample}.bam'),
        os.path.join(OD_BAM,'{sample}.bam.bai')
    output:
        temp(os.path.join(OD_METRICS,'{sample}.ensembl.RNAmetrics.txt'))
    log:
        os.path.join(OD_LOG,'{sample}.picard.log')
    params:
        str = PICARD_STRAND
    threads: 1
    resources:
        mem_mb = 20000
    singularity:
        config['PICARD_IMAGE']
    shell:
        """
        export _JAVA_OPTIONS="-Xmx10g" && \
        /usr/local/bin/picard CollectRnaSeqMetrics \
            --REF_FLAT {REFFLAT_ENS} \
            --RIBOSOMAL_INTERVALS {RIBO_INTERVALS} \
            --STRAND_SPECIFICITY {params.str} \
            -I {input[0]} \
            -O {output} 2> {log}
        """

# ------------------------------------------------------------------------------
rule picard_refseq:
    input:
        os.path.join(OD_BAM,'{sample}.bam'),
        os.path.join(OD_BAM,'{sample}.bam.bai')
    output:
        temp(os.path.join(OD_METRICS,'{sample}.refseq.RNAmetrics.txt'))
    log:
        os.path.join(OD_LOG,'{sample}.refseq.picard.log')
    params:
        str = PICARD_STRAND
    threads: 1
    resources:
        mem_mb = 20000
    singularity:
        config['PICARD_IMAGE']
    shell:
        """
        export _JAVA_OPTIONS="-Xmx10g" && \
        /usr/local/bin/picard CollectRnaSeqMetrics \
            --REF_FLAT {REFFLAT_REF} \
            --RIBOSOMAL_INTERVALS {RIBO_INTERVALS} \
            --STRAND_SPECIFICITY {params.str} \
            -I {input[0]} \
            -O {output} 2> {log}
        """
