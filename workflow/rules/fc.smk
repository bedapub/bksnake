"""
Calculate gene counts with feature counts from subreads
"""

# ------------------------------------------------------------------------------
# Variables

MIN_OVERLAP = config['fc_min_overlap']
FRAC_OVERLAP = config['fc_frac_overlap']

"""
NOTE: Introduction of new parameter in v2.0.2 (https://subread.sourceforge.net): 
  
  --countReadPairs

New parameter '--countReadPairs' is added to featureCounts to explicitly specify that read pairs will be counted, 
and the '-p' option in featureCounts now only specifies if the input reads are paired end 
(it also implied that counting of read pairs would be performed in previous versions).
"""
if config['library']['type'] == 'paired-end':
    P_PARAM = '-p --countReadPairs'
else:
    P_PARAM = ''


# ------------------------------------------------------------------------------
rule fc:
    input:
        bam = os.path.join(OD_UBAM,'{sample}_Aligned.out.bam'),
        gtf = os.path.join(OD_ANNO,'{db}.gtf.gz')
    output:
        os.path.join(OD_FC,'{sample}.{db}.cnt.gz'),
        os.path.join(OD_FC,'{sample}.{db}.cnt.summary')
    log:
        os.path.join(OD_LOG,'{sample}.fc_{db}.log')
    params:
        cnt = os.path.join(OD_FC,'{sample}.{db}.cnt')
    threads: 4
    resources:
        mem_mb = 30000
    singularity:
        config['SUBREAD_IMAGE']
    shell:
        """
        featureCounts \
            -t exon -g gene_id {P_PARAM} \
            -Q 10 -B -C --minOverlap {MIN_OVERLAP} \
             --fracOverlap {FRAC_OVERLAP} \
            -s {FC_STRAND} -a {input.gtf} \
            -T {threads} -o {params.cnt} {input.bam} 2> {log} \
        && gzip {params.cnt}
        """