RIBIOS_IMAGE = config['RIBIOS_IMAGE']
CBIOS_IMAGE = config['CBIOS_IMAGE']

rule gct_ref:
    input:
        expand(os.path.join(OD_FC,'{sample}.refseq.cnt.gz'), sample=sample_ids)
    output:
        {COUNT_GCT_REF}
    log:
        os.path.join(OD_LOG,'gct_ref.log')
    params:
        annot = {ANNOT_REF},
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
            --annotation {params.annot} \
            {input} {output} 2> {log}
        """

rule gct_ens:
    input:
        expand(os.path.join(OD_FC,'{sample}.ensembl.cnt.gz'), sample=sample_ids)
    output:
        {COUNT_GCT_ENS}
    log:
        os.path.join(OD_LOG,'gct_ens.log')
    params:
        annot = {ANNOT_ENS},
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
            --annotation {params.annot} \
            {input} {output} 2> {log}
        """

rule tpm_ref:
    input:
        {COUNT_GCT_REF}
    output:
        {TPM_GCT_REF},
        {LOG2TPM_GCT_REF}
    log:
        os.path.join(OD_LOG,'tpm_ref.log')
    params:
        length = {LEN_REF},
        col = 4
    threads: 1
    resources:
        mem_mb = 10000
    singularity:
        CBIOS_IMAGE
    shell:
        """
        count2tpm \
            -g {input} \
            -l {params.length} \
            -tpm -digits 3 \
            -col {params.col} > {output[0]} 2> {log}
        count2tpm \
            -log2 \
            -g {input} \
            -l {params.length} \
            -tpm -digits 6 \
            -col {params.col} > {output[1]} 2>> {log}
        """

rule tpm_ens:
    input:
        {COUNT_GCT_ENS}
    output:
        {TPM_GCT_ENS},
        {LOG2TPM_GCT_ENS}
    log:
        os.path.join(OD_LOG,'tpm_ens.log')
    params:
        length = {LEN_ENS},
        col = 4
    threads: 1
    resources:
        mem_mb = 10000
    singularity:
        CBIOS_IMAGE
    shell:
        """
        count2tpm \
            -g {input} \
            -l {params.length} \
            -tpm -digits 3 \
            -col {params.col} > {output[0]} 2> {log}
        count2tpm \
            -log2 \
            -g {input} \
            -l {params.length} \
            -tpm -digits 6 \
            -col {params.col} > {output[1]} 2>> {log}
         """

# Convert species genes to human genes for gct matrix: count
rule collapse:
    input:
        ref = {COUNT_GCT_REF},
        ens = {COUNT_GCT_ENS},
    output:
        ref = {COLLAPSED_GCT_REF}, 
        ens = {COLLAPSED_GCT_ENS}, 
    params:
        chipFile_ref = {CHIP_FILE_REF},
        chipFile_ens = {CHIP_FILE_ENS}
    resources:
        mem_mb = 10000
    singularity:
        RIBIOS_IMAGE
    shell:
        """
        collapseExprsMatByChip.Rscript -useChipfileAnno \
            -infile {input.ref} -outfile {output.ref} -chipfile {params.chipFile_ref}
        collapseExprsMatByChip.Rscript -useChipfileAnno \
            -infile {input.ens} -outfile {output.ens} -chipfile {params.chipFile_ens}
        """

# Convert species genes to human genes for gct matrix: tpm
rule collapse_tpm:
    input:
        ref = {TPM_GCT_REF},
        ens = {TPM_GCT_ENS},
    output:
        ref = {TPM_COLLAPSED_GCT_REF},
        ens = {TPM_COLLAPSED_GCT_ENS},
    params:
        chipFile_ref = {CHIP_FILE_REF},
        chipFile_ens = {CHIP_FILE_ENS}
    resources:
        mem_mb = 10000
    singularity:
        RIBIOS_IMAGE
    shell:
        """
        collapseExprsMatByChip.Rscript -useChipfileAnno \
            -infile {input.ref} -outfile {output.ref} -chipfile {params.chipFile_ref}
        collapseExprsMatByChip.Rscript -useChipfileAnno \
            -infile {input.ens} -outfile {output.ens} -chipfile {params.chipFile_ens}
        """

# Traditional QC files (bioQC, PCA)
rule qc:
    input:
        pheno = {PHENODATA},
        tpm_ref = {TPM_GCT_REF},
        collapsed = {COLLAPSED_GCT_REF}, 
        tpm_ens = {TPM_GCT_ENS},
        log2_ref = {LOG2TPM_GCT_REF},
        log2_ens = {LOG2TPM_GCT_ENS}
    output:
        cls = os.path.join(OD_ANNO,'phenoData.cls'),
        pca = os.path.join(OD_QC,'refseq_log2tpm_pca.pdf'),
        tab = os.path.join(OD_QC,'refseq_log2tpm_pca.txt'),
        txt = os.path.join(OD_QC,'bioQC.txt'),
        txt2 = os.path.join(OD_QC,'bioQC_thr2.txt'),
        pdf = os.path.join(OD_QC,'bioQC.pdf'),
    log:
        os.path.join(OD_LOG,'qc.log')
    params:
        gmt = 'resources/exp.tissuemark.bioqc.roche.symbols.gmt'
    resources:
        mem_mb = 10000
    singularity:
        RIBIOS_IMAGE
    shell:
        """
        R -e \"df=ribiosIO::readTable(\'{input.pheno}\'); ribiosIO::write_cls(factor(df\$GROUP),\'{output.cls}\')\" > {log} 
        (plotPCA.Rscript -infile {input.log2_ref} -outfile {output.pca} -outTable {output.tab} -cls {output.cls}) >> {log}
        (bioqc.Rscript -gmt {params.gmt} -featuretype GeneSymbol -infile {input.collapsed} -outfile {output.txt}) >> {log}
        (biosHeatmap.Rscript -infile {output.txt} -outfile {output.pdf} -scale none -colors blackred \
            -naColor lightgray -symbreaks auto -dendrogram both -dist euclidean -hclust ward.D2 \
            -xlab -ylab -cexRow -cexCol -colorKeyTitle -width -height -margins -zlimLo -zlimHi -main \'BioQC\') >> {log}
        (bioqc.Rscript -threshold 2.0 -gmt {params.gmt} -featuretype GeneSymbol -infile {input.collapsed} -outfile {output.txt2}) >> {log}
        """
