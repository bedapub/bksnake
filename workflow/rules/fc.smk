# TO DO: modularize better, merge ref and ens

MIN_OVERLAP = config['fc_min_overlap']
FRAC_OVERLAP = config['fc_frac_overlap']
if config['library']['type'] == 'paired-end':
    P_PARAM = '-p'
else:
    P_PARAM = ''


# ------------------------------------------------------------------------------
rule fc_ref:
    input:
        os.path.join(OD_UBAM,'{sample}_Aligned.out.bam')
    output:
        os.path.join(OD_FC,'{sample}.refseq.cnt.gz'),
        os.path.join(OD_FC,'{sample}.refseq.cnt.summary')
    log:
        os.path.join(OD_LOG,'{sample}.fc_ref.log')
    params:
        cnt = os.path.join(OD_FC,'{sample}.refseq.cnt')
    threads: 1
    resources:
        mem_mb = 30000
    singularity:
        config['SUBREAD_IMAGE']
    shell:
        'featureCounts '
        '    -t exon -g gene_id {P_PARAM}'
        '    -Q 10 -B -C --minOverlap {MIN_OVERLAP}'
        '     --fracOverlap {FRAC_OVERLAP}'
        '    -s {FC_STRAND} -a {GTF_REF}'
        '    -T {threads} -o {params.cnt} {input} 2> {log} &&'
        'gzip {params.cnt}'


# ------------------------------------------------------------------------------
rule fc_ens:
    input:
        os.path.join(OD_UBAM,'{sample}_Aligned.out.bam')
    output:
        os.path.join(OD_FC,'{sample}.ensembl.cnt.gz'),
        os.path.join(OD_FC,'{sample}.ensembl.cnt.summary')
    log:
        os.path.join(OD_LOG,'{sample}.fc_ens.log')
    params:
        cnt = os.path.join(OD_FC,'{sample}.ensembl.cnt')
    threads: 1
    resources:
        mem_mb = 30000
    singularity:
        config['SUBREAD_IMAGE']
    shell:
        'featureCounts '
        '    -t exon -g gene_id {P_PARAM}'
        '    -Q 10 -B -C --minOverlap {MIN_OVERLAP}'
        '    --fracOverlap {FRAC_OVERLAP}'
        '    -s {FC_STRAND} -a {GTF_ENS}'
        '    -T {threads} -o {params.cnt} {input} 2> {log} &&'
        'gzip {params.cnt}'
