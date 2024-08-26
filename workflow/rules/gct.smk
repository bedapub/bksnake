# ------------------------------------------------------------------------------
# Variables 
RIBIOSSCRIPTS_IMAGE = config['RIBIOSSCRIPTS_IMAGE']
NGSTOOLS_IMAGE = config['NGSTOOLS_IMAGE']


# ------------------------------------------------------------------------------
"""
Create GCT file with count for RefSeq/Ensembl gene annotations
"""
rule gct:
    input:
        files = lambda wildcards: expand(rules.fc.output.cnt, sample=sample_ids, db=wildcards.db),
        annot = rules.annotations.output.annot,
    output:
        gct = os.path.join(OD_GCT, '{db}_count.gct'),
    log:
        os.path.join(OD_LOG, 'gct_{db}.log')
    group: 'gct'
    params:
        order = ','.join(sample_ids)
    threads: 1
    resources:
        mem_mb = 10000
    shell:
        """
        python workflow/scripts/fc2gct.py \
            --remove-suffix-string _Aligned.out.bam \
            --basename \
            --order {params.order} \
            --annotation {input.annot} \
            {input.files} {output} 2> {log}
        """

# ------------------------------------------------------------------------------
"""
Calculate normalized counts, TPM from RefSeq/Ensembl counts GCT file
"""
rule tpm:
    input:
        gct = rules.gct.output.gct,
        len = rules.annotations.output.len,
    output:
        tpm = os.path.join(OD_GCT, '{db}_tpm.gct'),
        log = os.path.join(OD_GCT, '{db}_log2tpm.gct')
    log:
        os.path.join(OD_LOG, 'tpm_{db}.log')
    group: 'gct'
    params:
        col = GENE_LENGTH_COLUMN_INDEX,
        len = os.path.join(OD_GCT, '{db}.geneLength')
    threads: 1
    resources:
        mem_mb = 10000
    singularity:
        NGSTOOLS_IMAGE
    shell:
        """
        gunzip -c {input.len} > {params.len}
        count2tpm \
            -g {input.gct} \
            -l {params.len} \
            -tpm -digits 3 \
            -col {params.col} > {output.tpm} 2> {log}
        count2tpm \
            -log2 \
            -g {input.gct} \
            -l {params.len} \
            -tpm -digits 6 \
            -col {params.col} > {output.log} 2>> {log}
        rm -f {params.len}
        """

# ------------------------------------------------------------------------------
"""
Convert species genes to human genes for gct matrix: count and tpm
Uses Chip definition files

Note: the gencode GCT file contain currenlty a dot in the gene id, e.g. "ENSG00000000003.16"
temporarily th number will be removed for the "collapse" rule
        #chip = os.path.join('resources', '{db}.chip')

"""
rule collapse:
    input:
        tpm = rules.tpm.output.tpm,
        cnt = rules.gct.output.gct,
    output:
        tpm = os.path.join(OD_GCT, '{db}_tpm_collapsed.gct'),
        cnt = os.path.join(OD_GCT, '{db}_count_collapsed.gct')
    log:
        os.path.join(OD_LOG, 'collapse_{db}.log')
    group: 'gct'
    params:
        chip = lambda wildcards: os.path.join('resources', f"{wildcards.db.split('_')[0] if '_' in wildcards.db else wildcards.db}.chip")
    resources:
        mem_mb = 10000
    singularity:
        RIBIOSSCRIPTS_IMAGE
    shell:
        """
        collapseExprsMatByChip.Rscript -useChipfileAnno \
            -infile {input.tpm} -outfile {output.tpm} -chipfile {params.chip}
        collapseExprsMatByChip.Rscript -useChipfileAnno \
            -infile {input.cnt} -outfile {output.cnt} -chipfile {params.chip}
        """
        
# ------------------------------------------------------------------------------
"""
Create traditional QC files (bioQC, PCA plots) based on RefSeq or Ensembl.
     labels={
              "model": "{model}",
              "figure": "some plot"
          }
"""
rule qc:
    input:
        pheno = {PHENODATA},
        log2 = rules.tpm.output.log,
        collapsed = rules.collapse.output.tpm,
    output:
        cls = temp(os.path.join(OD_ANNO, '{db}_phenoData.cls')),
        pca = report(os.path.join(OD_QC, '{db}_log2tpm_pca.png'), 
            category='Quality Control', subcategory='{db}', 
            caption='../report/fig1_pca.rst',
            labels={
              'annotation': '{db}',
              'figure': 'Plot of the principal component analysis'
              }),
        tab = os.path.join(OD_QC, '{db}_log2tpm_pca.txt'),
        txt = os.path.join(OD_QC, '{db}_bioQC.txt'),
        thr = os.path.join(OD_QC, '{db}_bioQC_thr2.txt'),
        bioqc = report(os.path.join(OD_QC, '{db}_bioQC.png'), 
            category='Quality Control', subcategory='{db}', 
            caption='../report/fig2_bioqc.rst',
            labels={
              'annotation': '{db}',
              'figure': 'BioQC tissue heterogeneity plot'
              }),
    log:
        os.path.join(OD_LOG, '{db}_qc.log')
    params:
    resources:
        mem_mb = 10000
    singularity:
        RIBIOSSCRIPTS_IMAGE
    shell:
        """
        R -e \"df=ribiosIO::readTable(\'{input.pheno}\'); ribiosIO::write_cls(factor(df\$GROUP),\'{output.cls}\')\" > {log}
        
        (plotPCA.Rscript -infile {input.log2} -outfile {output.pca} -outTable {output.tab} -cls {output.cls}) >> {log}
        
        (bioqc.Rscript -featuretype GeneSymbol -infile {input.collapsed} -outfile {output.txt}) >> {log}
        
        (biosHeatmap.Rscript -infile {output.txt} -outfile {output.bioqc} -scale none -colors blackred \
            -naColor lightgray -symbreaks auto -dendrogram both -dist euclidean -hclust ward.D2 \
            -xlab -ylab -cexRow -cexCol -colorKeyTitle -width -height -margins -zlimLo -zlimHi -main \'BioQC\') >> {log}
        
        (bioqc.Rscript -threshold 2.0 -featuretype GeneSymbol -infile {input.collapsed} -outfile {output.thr}) >> {log}
        """
