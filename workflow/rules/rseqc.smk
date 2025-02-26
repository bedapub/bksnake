# -------------------------------------------------------------
rule strandness:
    input:
        bam = rules.sortbam_star.output,
        bai = rules.indexbam.output,
	bed = rules.annotations.output.bed,
    output:
        txt = temp(os.path.join(OD_METRICS, '{sample}.{db}.strandness.txt')),
        tmp1 = temp(os.path.join(OD_METRICS, '{sample}.{db}.strandness.tmp1')),
        tmp2 = temp(os.path.join(OD_METRICS, '{sample}.{db}.strandness.tmp2')),
        bed = temp(os.path.join(OD_METRICS, '{sample}.{db}.bed')),
    log:
        os.path.join(OD_LOG, '{sample}.{db}.strandness.log')
    params:
        strandness = config['strandness_mode'],
    threads: 1
    resources:
        mem_mb = 10000
    singularity:
        config['RSEQC_IMAGE']
    shell:
        """
        gunzip -c {input.bed} > {output.bed} \
        && infer_experiment.py -r {output.bed} -i {input.bam} > {output.tmp1} \
        && workflow/scripts/strandness.sh {output.tmp1} {params.strandness} > {output.tmp2} \
	&& cat {output.tmp1} {output.tmp2} > {output.txt}
        """        

# -------------------------------------------------------------
rule gtfToGenePred:
    input:
        rules.annotations.output.gtf,
    output:
        temp(os.path.join(OD_ANNO, '{db}.gtfToGenePred')),
    log:
        os.path.join(OD_LOG, '{db}.gtfToGenePred.log')
    group: 'bedfile'
    threads: 1
    resources:
        mem_mb = 10000
    singularity:
        config['GTFTOGENEPRED_IMAGE']
    shell:
        """
        gzip -cd {input} | grep -v 'unknown_transcript_' | gtfToGenePred -ignoreGroupsWithoutExons /dev/stdin {output}
        """        

# -------------------------------------------------------------
rule genePredToBed:
    input:
        rules.gtfToGenePred.output,
    output:
        temp(os.path.join(OD_ANNO, '{db}.genePredToBed.bed')),
    log:
        os.path.join(OD_LOG, '{db}.genePredToBed.log')
    group: 'bedfile'
    threads: 1
    resources:
        mem_mb = 10000
    singularity:
        config['GENEPREDTOBED_IMAGE']
    shell:
        """
        genePredToBed {input} {output}
        """        

# -------------------------------------------------------------
rule junction_annotation:
    input:
        bam = rules.sortbam_star.output,
        bai = rules.indexbam.output,
        bed = rules.genePredToBed.output,
    output:
        os.path.join(OD_METRICS, '{sample}.{db}.junction_annotation.log'),
        os.path.join(OD_METRICS, '{sample}.{db}.splice_events.pdf'),
        os.path.join(OD_METRICS, '{sample}.{db}.splice_junction.pdf'),
        os.path.join(OD_METRICS, '{sample}.{db}.junction_plot.r'),
        os.path.join(OD_METRICS, '{sample}.{db}.junction.xls'),
        os.path.join(OD_METRICS, '{sample}.{db}.junction.bed'),
        os.path.join(OD_METRICS, '{sample}.{db}.junction.Interact.bed'),
    log:
        os.path.join(OD_LOG, '{sample}.{db}.junction_annotation.log')
    params:
        prefix = os.path.join(OD_METRICS, '{sample}.{db}')
    threads: 1
    resources:
        mem_mb = 10000
    singularity:
        config['RSEQC_IMAGE']
    shell:
        """
        junction_annotation.py -i {input.bam} -r {input.bed} -o {params.prefix} 2> {output[0]}
        """        

# -------------------------------------------------------------
rule junction_saturation:
    input:
        bam = rules.sortbam_star.output,
        bai = rules.indexbam.output,
        bed = rules.genePredToBed.output,
    output:
        os.path.join(OD_METRICS, '{sample}.{db}.junctionSaturation_plot.r'),
        os.path.join(OD_METRICS, '{sample}.{db}.junctionSaturation_plot.pdf'),
    log:
        os.path.join(OD_LOG, '{sample}.{db}.junction_saturation.log')
    params:
        prefix = os.path.join(OD_METRICS, '{sample}.{db}')
    threads: 1
    resources:
        mem_mb = 10000
    singularity:
        config['RSEQC_IMAGE']
    shell:
        """
        junction_saturation.py -i {input.bam} -r {input.bed} -o {params.prefix} > {log}
        """
        
