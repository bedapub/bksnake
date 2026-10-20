# ------------------------------------------------------------------------------
# Picard mRNA metrics only for paired-end reads and for RefSeq/Ensembl annotations
#        ref = rules.genome.output.ugz,

rule picard:
    input:
        bam = rules.sortbam_star.output,
        bai = rules.indexbam.output,
        ref = rules.annotations.output.flat,
        ribo = rules.annotations.output.ribo,
        str = rules.strandness.output.txt,
    output:
        temp(os.path.join(OD_METRICS, '{sample}.{db}.RNAmetrics.txt'))
    log:
        os.path.join(OD_LOG, '{sample}.{db}.picard.log')
    params:
        records = 10000000, # default is 500000, only used by Picard when writing BAM files, no effect on runtime or memory here
    threads: config.get('picard_threads', 8)
    resources:
        mem_mb = 30000
    singularity:
        config['PICARD_IMAGE']
    shell:
        """
        str="$(grep 'picard=' {input.str} | sort | uniq | cut -d= -f2)"
        
        export _JAVA_OPTIONS="-Xmx25g" && \
        workflow/scripts/picard_rnaseq_metrics.sh {threads} \
            --REF_FLAT {input.ref} \
            --RIBOSOMAL_INTERVALS {input.ribo} \
            --STRAND_SPECIFICITY ${{str}} \
            --MAX_RECORDS_IN_RAM {params.records} \
            -I {input.bam} \
            -O {output} 2> {log}
        """
