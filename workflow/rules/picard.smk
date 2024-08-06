# -------------------------------------------------------------
rule strandedness:
    input:
        bam = os.path.join(OD_BAM, '{sample}.bam'),
        bai = os.path.join(OD_BAM, '{sample}.bam.bai')
    output:
        txt = temp(os.path.join(OD_METRICS, '{sample}.strandedness.txt')),
        bed = temp(os.path.join(OD_METRICS, '{sample}.bed')),
    log:
        os.path.join(OD_LOG, '{sample}.strandedness.log')
    params:
        bed = os.path.join(OD_ANNO, DBS[0]+'.bed.gz'),
    threads: 1
    resources:
        mem_mb = 10000
    singularity:
        config['RSEQC_IMAGE']
    shell:
        """
        gunzip -c {params.bed} > {output.bed}
        infer_experiment.py -r {output.bed} -i {input.bam} > {output.txt}
        python workflow/scripts/strandedness.py {output.txt} >> {output.txt}
        """        


# ------------------------------------------------------------------------------
# Picard mRNA metrics only for paired-end reads and for RefSeq/Ensembl annotations
rule picard:
    input:
        bam = os.path.join(OD_BAM, '{sample}.bam'),
        bai = os.path.join(OD_BAM, '{sample}.bam.bai'),
        ref = os.path.join(OD_ANNO, '{db}.refFlat.gz'),
        ribo = os.path.join(OD_ANNO, 'genome.rRNA_intervals'),
        str = os.path.join(OD_METRICS, '{sample}.strandedness.txt'),
    output:
        temp(os.path.join(OD_METRICS, '{sample}.{db}.RNAmetrics.txt'))
    log:
        os.path.join(OD_LOG, '{sample}.{db}.picard.log')
    params:
        records = 10000000, # default is 500000
    threads: 1
    resources:
        mem_mb = 30000
    singularity:
        config['PICARD_IMAGE']
    shell:
        """
        str="$(grep 'picard=' {input.str} | cut -d= -f2)"
        
        export _JAVA_OPTIONS="-Xmx25g" && \
        /usr/local/bin/picard CollectRnaSeqMetrics \
            --REF_FLAT {input.ref} \
            --RIBOSOMAL_INTERVALS {input.ribo} \
            --STRAND_SPECIFICITY ${{str}} \
            --MAX_RECORDS_IN_RAM {params.records} \
            -I {input.bam} \
            -O {output} 2> {log}
        """
