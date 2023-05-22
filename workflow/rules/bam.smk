SAMTOOLS_IMAGE = config['SAMTOOLS_IMAGE']

# -------------------------------------------------------------
rule indexbam:
    input:
        os.path.join(OD_BAM, '{sample}.bam'),
    output:
        (os.path.join(OD_BAM, '{sample}.bam.bai'))
    threads: 2
    resources:
        mem_mb = 20000
    singularity:
         SAMTOOLS_IMAGE
    shell:
        'samtools index -@ {threads} {input}'

# -------------------------------------------------------------
rule flagstat:
    input:
        os.path.join(OD_BAM, '{sample}.bam')
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
        os.path.join(OD_BAM, '{sample}.bam')
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
        os.path.join(OD_BAM, '{sample}.bam')
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
        os.path.join(OD_BAM, '{sample}.bam')
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
        os.path.join(OD, 'rmdup/{sample}.bam')
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
        bam = os.path.join(OD_BAM, '{sample}.bam'),
        bai = os.path.join(OD_BAM, '{sample}.bam.bai'),
        genome_fa = os.path.join(OD_ANNO, 'genome.fa.gz'),
        genome_fa_fai = os.path.join(OD_ANNO, 'genome.fa.gz.fai'),
        genome_fa_gzi = os.path.join(OD_ANNO, 'genome.fa.gz.gzi')
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
        os.path.join(OD_CRAM, '{sample}.cram')
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
# only data from flagstat
rule mapping_stats:
    input:
        expand(os.path.join(OD_STATS, '{sample}.bam.flagstat'), sample=sample_ids)
    output:
        os.path.join(OD_STATS, 'flagstat_mapping_stats.txt')
    threads: 1
    resources:
        mem_mb = 1000
    run:
        with open(output[0], 'w') as f:
            f.write('ID\tGROUP\tTOTAL_READS\tMAPPED_READS\tMAPPED_IN_PERC\tUNMAPPED_READS\tUNMAPPED_IN_PERC\n')
#            f.write('SampleName\tID\tGROUP\tTOTAL_READS\tMAPPED_READS\tMAPPED_IN_PERC\tUNMAPPED_READS\tUNMAPPED_IN_PERC\n')
            for ind, sid in enumerate(sample_ids):
#                f.write(sample_names[ind] + '\t')
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
# bigwig files for jbrowse
# Documentation for bamCoverage
# https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html

if 'jbrowse' in config:
    binSize = config['jbrowse']['bamCoverage_binSize']
else:
    binSize = 3

if config['generate_bw_files'] == True:
    rule bw:
        input:
            os.path.join(OD_BAM, '{sample}.bam'),
            os.path.join(OD_BAM, '{sample}.bam.bai')
        output:
            os.path.join(OD_BW, '{sample}.bw'),
            temp(os.path.join(OD_BW, '{sample}.bw.done'))
        log:
            os.path.join(OD_LOG, '{sample}.bw.log')
        params:
            prefix = os.path.join(OD_BW, '{sample}'),
            binSize = binSize
        threads: 8
        resources:
            mem_mb = 20000
        singularity:
            config['DEEPTOOLS_IMAGE']
        shell:
            """
            (bamCoverage --ignoreDuplicates --binSize {params.binSize} \
                -p {threads} -b {input[0]} -o {output[0]}) 2> {log} && \
            touch {output[1]}
            """
else:
    rule touch_bw:
        output:
            temp(os.path.join(OD_BW,'{sample}.bw.done'))
        shell:
            """
            touch {output}
            """                        
