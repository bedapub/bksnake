SAMTOOLS_IMAGE = config['SAMTOOLS_IMAGE']

# -------------------------------------------------------------
rule indexbam:
    input:
        rules.sortbam_star.output,
    output:
        os.path.join(OD_BAM, '{sample}.bam.bai') if config['keep_bam_files'] else temp(os.path.join(OD_BAM, '{sample}.bam.bai')),
    threads: 2
    resources:
        mem_mb=20000
    singularity:
        SAMTOOLS_IMAGE
    shell:
        'samtools index -@ {threads} {input}'
        
# -------------------------------------------------------------
rule flagstat:
    input:
        rules.sortbam_star.output,
    output:
        temp(os.path.join(OD_STATS, '{sample}.bam.flagstat'))
    threads: 2
    resources:
        mem_mb = 20000
    singularity:
         SAMTOOLS_IMAGE
    shell:
        'samtools flagstat -@ {threads} {input} > {output}'

# -------------------------------------------------------------
rule samstats:
    input:
        rules.sortbam_star.output,
    output:
        temp(os.path.join(OD_STATS, '{sample}.bam.stats'))
    threads: 2
    resources:
        mem_mb = 20000
    singularity:
         SAMTOOLS_IMAGE
    shell:
        'samtools stats -@ {threads} {input} > {output}'

# -------------------------------------------------------------
rule bamstats:
    input:
        rules.sortbam_star.output,        
    output:
        temp(os.path.join(OD_STATS, '{sample}.bam.stats2'))
    threads: 1
    resources:
        mem_mb = 20000
    singularity:
        config['BAMTOOLS_IMAGE']
    shell:
        'bamtools stats -in {input} > {output}'

# -------------------------------------------------------------
rule rmdup:
    input:
        rules.sortbam_star.output,
    output:
        temp(os.path.join(OD, 'rmdup/{sample}.bam'))
    threads: 1
    resources:
        mem_mb = 20000
    singularity:
        SAMTOOLS_IMAGE
    shell:
        'samtools rmdup {input} {output}'

# -------------------------------------------------------------
rule indexrmdup:
    input:
        rules.rmdup.output,
    output:
        temp(os.path.join(OD, 'rmdup/{sample}.bam.bai'))
    threads: 2
    resources:
        mem_mb = 20000
    singularity:
        SAMTOOLS_IMAGE
    shell:
        'samtools index -@ {threads} {input}'

# -------------------------------------------------------------
rule cram:
    input:
        bam = rules.sortbam_star.output,
        bai = rules.indexbam.output,
        genome_fa = rules.genome.output.gz,
        fai = rules.genome.output.fai,
        gzi = rules.genome.output.gzi,
    output:
        os.path.join(OD_CRAM, '{sample}.cram')
    log:
        os.path.join(OD_LOG, '{sample}.cram.log')
    threads: 4
    resources:
        mem_mb = 20000
    singularity:
         SAMTOOLS_IMAGE
    shell:
        'samtools view -@ {threads} -C -T {input.genome_fa} {input.bam} -O CRAM -o {output} 2> {log}'

# -------------------------------------------------------------
rule cramindex:
    input:
        rules.cram.output,
    output:
        os.path.join(OD_CRAM, '{sample}.cram.crai')
    log:
        os.path.join(OD_LOG, '{sample}.cramindex.log')
    threads: 2
    resources:
        mem_mb = 10000
    singularity:
         SAMTOOLS_IMAGE
    shell:
        'samtools index -@ {threads} {input} 2> {log}'

# -------------------------------------------------------------
# Use only data from flagstat
rule mapping_stats:
    input:
        expand(rules.flagstat.output, sample=sample_ids),
    output:
        os.path.join(OD_STATS, 'flagstat_mapping_stats.txt')
    threads: 1
    resources:
        mem_mb = 1000
    run:
        with open(output[0], 'w') as f:
            f.write('ID\tGROUP\tTOTAL_READS\tMAPPED_READS\tMAPPED_IN_PERC\tUNMAPPED_READS\tUNMAPPED_IN_PERC\n')
            for ind, sid in enumerate(sample_ids):
                f.write(str(sid) + '\t')
                f.write(sample_groups[ind] + '\t')
                with open(input[ind], 'r') as fsf:
                    total_reads = fsf.readline().split(' ')[0]
                    for _ in range(3):
                        next(fsf)
                    mapped_reads = fsf.readline().split(' ')[0]
                    mapping_rate = int(mapped_reads)/int(total_reads)
                    unmapped_reads = int(total_reads) - int(mapped_reads)
                    unmapping_rate = int(unmapped_reads)/int(total_reads)
                    f.write('{}\t{}\t{:.3f}\t{}\t{:.3f}\n'.format(
                        total_reads,
                        mapped_reads,
                        mapping_rate,
                        unmapped_reads,
                        unmapping_rate))
                        
# -------------------------------------------------------------
# BigWig files for JBrowse genome viewer
# Documentation for bamCoverage:
#   https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
rule bw:
    input:
        rules.sortbam_star.output,
        rules.indexbam.output,
    output:
        os.path.join(OD_BW, '{sample}.bw'),
    log:
        os.path.join(OD_LOG, '{sample}.bw.log')
    params:
        prefix = os.path.join(OD_BW, '{sample}'),
        binSize = config['jbrowse']['bamCoverage_binSize'] if 'jbrowse' in config and 'bamCoverage_binSize' in config['jbrowse'] else 3,
    threads: 8
    resources:
        mem_mb = 20000
    singularity:
        config['DEEPTOOLS_IMAGE']
    shell:
        """
        (bamCoverage --ignoreDuplicates --binSize {params.binSize} \
            -p {threads} -b {input[0]} -o {output}) 2> {log}
        """
