# ------------------------------------------------------------------------------
# Picard mRNA metrics only for paired-end reads and for RefSeq/Ensembl annotations
rule picard:
    input:
        bam = os.path.join(OD_BAM, '{sample}.bam'),
        bai = os.path.join(OD_BAM, '{sample}.bam.bai'),
        ref = os.path.join(OD_ANNO, '{db}.refFlat.gz'),
        ribo = os.path.join(OD_ANNO, 'genome.rRNA_intervals')
    output:
        temp(os.path.join(OD_METRICS, '{sample}.{db}.RNAmetrics.txt'))
    log:
        os.path.join(OD_LOG, '{sample}.{db}.picard.log')
    params:
        str = PICARD_STRAND, 
        records = 1000000 # default is 500000 
    threads: 1
    resources:
        mem_mb = 30000
    singularity:
        config['PICARD_IMAGE']
    shell:
        """
        export _JAVA_OPTIONS="-Xmx25g" && \
        /usr/local/bin/picard CollectRnaSeqMetrics \
            --REF_FLAT {input.ref} \
            --RIBOSOMAL_INTERVALS {input.ribo} \
            --STRAND_SPECIFICITY {params.str} \
            --MAX_RECORDS_IN_RAM {params.records} \
            -I {input.bam} \
            -O {output} 2> {log}
        """
