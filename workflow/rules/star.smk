"""--------------------------------------------------------------------
Function to create a data table with mapping statistics from STAR
"""
def star_mapping_stats (outfile, phenoFile, indir):
    """Function to create a data table with mapping statistics from STAR
    Parameters
    ----------
    outfile : str
        Name of output file.
    phenoFile : str
        Name of file containing sample metadata
    indir : str
        Name of directory containing STAR _Log.final.out file(s)
    """
    pheno = pd.read_csv(phenoFile, header = 0, sep = '\t', usecols=['ID','GROUP'],
            dtype={'ID':str,'GROUP':str},
            keep_default_na=False, na_values=[''])
    # Read in all files with star output
    i = 0
    for sample in list(pheno['ID']):
        f = os.path.join(indir, sample+'_Log.final.out')
        df = pd.read_csv(f, header = None, sep = '\t')
        df.columns = ['0', sample]
        if i == 0:
            tab = df
        else:
            tab = pd.concat([tab,df[sample]], axis = 1)
        i = 1
    del df

    # Re-format data table and merge with annotation data table
    tab['0'] = tab['0'].str.replace(r'\|', '', regex=True)
    tab['0'] = tab['0'].str.strip()
    df = tab.transpose()
    df.reset_index(inplace=True)
    df.columns = df.iloc[0]
    df.drop(df.index[0], inplace = True)
    df.rename(columns={"0": "ID"}, inplace = True)

    # Calculated mapping rates
    tab = pd.merge(df, pheno, on=['ID'])
    tab['TOTAL_READS'] = pd.to_numeric(tab['Number of input reads'])
    tab['UNMAPPED_READS'] = pd.to_numeric(tab['Number of reads unmapped: too many mismatches']) + \
                            pd.to_numeric(tab['Number of reads unmapped: too short']) + \
                            pd.to_numeric(tab['Number of reads unmapped: other']) + \
                            pd.to_numeric(tab['Number of chimeric reads'])
    tab['MAPPED_READS'] = tab['TOTAL_READS'] - tab['UNMAPPED_READS']
    tab['MAPPED_IN_PERC'] = tab['MAPPED_READS'] / tab['TOTAL_READS']
    tab['UNMAPPED_IN_PERC'] = tab['UNMAPPED_READS'] / tab['TOTAL_READS']
    df = tab[['ID','GROUP','TOTAL_READS','MAPPED_READS','MAPPED_IN_PERC','UNMAPPED_READS','UNMAPPED_IN_PERC']]

    # Write table to output file
    df.to_csv(outfile, sep = '\t', index = False)


# ------------------------------------------------------------------------------
if config['generate_unmapped'] == True:
    unmapped = '--outReadsUnmapped Fastx'
else:
    unmapped = ''

if config['cutadapt']['run'] == True:
    star_input_dir = OD_CUTADAPT
else:
    star_input_dir = OD_FASTQ


# ------------------------------------------------------------------------------
# Try calculate memory for star per thread: we want in total 120Gb
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html
def get_mem_mb(wildcards, threads):
    return 2000 * threads # 20G for each threads
    #return 120 * 1024 / threads # 120G / threads -->FAILS


# ------------------------------------------------------------------------------
# paired-end read mapping with STAR
if config['library']['type'] == 'paired-end':
    rule star:
        input:
            fq1 = os.path.join(star_input_dir, '{sample}_1.fastq.gz'),
            fq2 = os.path.join(star_input_dir, '{sample}_2.fastq.gz'),
            cutadapt = os.path.join(OD_CUTADAPT, '{sample}.report.txt'),
            done1 = os.path.join(OD_FASTQ, '{sample}_1.fastq.gz.done'),
            done2 = os.path.join(OD_FASTQ, '{sample}_2.fastq.gz.done'),
        output:
            temp(directory(os.path.join(OD_UBAM, '{sample}__STARtmp'))),
            temp(directory(os.path.join(OD_UBAM, '{sample}__STARgenome'))),
            temp(os.path.join(OD_LOG, '{sample}_Log.final.out')),
            temp(os.path.join(OD_UBAM, '{sample}_Aligned.out.bam')),
            temp(os.path.join(OD_UBAM, '{sample}_Log.progress.out')),
            temp(os.path.join(OD_UBAM, '{sample}_SJ.out.tab')),
            temp(os.path.join(OD_UBAM, '{sample}_Unmapped.out.mate1')),
            temp(os.path.join(OD_UBAM, '{sample}_Unmapped.out.mate2'))
        log:
            f0 = os.path.join(OD_LOG, '{sample}.star.log'),
            f1 = os.path.join(OD_LOG, '{sample}_Log.out')
        retries: 3
        threads: 12 # config['star_threads']
        resources:
            mem_mb = get_mem_mb
            #mem_mb = 120 * 1024
        params:
            sam_mapq_unique = config['star_sam_mapq_unique'],
            prefix = join(OD_UBAM, '{sample}_'),
            rg = '{sample}'
        singularity:
            config['STAR_IMAGE']
        shell:
            """
            STAR --runThreadN {threads} \
                --genomeDir {STAR_DIR} \
                --outFileNamePrefix {params.prefix} \
                --outTmpDir {output[0]} \
                --outTmpKeep All \
                --sjdbGTFfile {GTF_FOR_STAR_MAPPING} \
                --outSAMmapqUnique {params.sam_mapq_unique} \
                --outSAMattributes All \
                --outSAMtype BAM Unsorted \
                --readFilesCommand zcat \
                --outReadsUnmapped Fastx \
                --outSAMattrRGline ID:{params.rg} SM:{params.rg} \
                --readFilesIn {input.fq1} {input.fq2} > {log.f0}
            mv --force {params.prefix}Log.out {log.f1}
            mv --force {params.prefix}Log.final.out {output[2]}
            """
else:
# ------------------------------------------------------------------------------
# single-end read mapping with STAR
    rule star:
        input:
            fq1 = os.path.join(star_input_dir, '{sample}.fastq.gz'),
            done1 = os.path.join(OD_FASTQ, '{sample}.fastq.gz.done'),
        output:
            temp(directory(os.path.join(OD_UBAM, '{sample}__STARtmp'))),
            temp(directory(os.path.join(OD_UBAM, '{sample}__STARgenome'))),
            temp(os.path.join(OD_LOG, '{sample}_Log.final.out')),
            temp(os.path.join(OD_UBAM, '{sample}_Aligned.out.bam')),
            temp(os.path.join(OD_UBAM, '{sample}_Log.progress.out')),
            temp(os.path.join(OD_UBAM, '{sample}_SJ.out.tab')),
            temp(os.path.join(OD_UBAM, '{sample}_Unmapped.out.mate1'))
        log:
            f0 = os.path.join(OD_LOG, '{sample}.star.log'),
            f1 = os.path.join(OD_LOG, '{sample}_Log.out')
        retries: 3
        threads: config['star_threads']
        resources:
            mem_mb = 120 * 1024
        params:
            sam_mapq_unique = config['star_sam_mapq_unique'],
            prefix = join(OD_UBAM, '{sample}_'),
	    rg = '{sample}'
        singularity:
            config['STAR_IMAGE']
        shell:
            """
            STAR --runThreadN {threads} \
                --genomeDir {STAR_DIR} \
                --outFileNamePrefix {params.prefix} \
                --outTmpDir {output[0]} \
                --outTmpKeep All \
                --sjdbGTFfile {GTF_FOR_STAR_MAPPING} \
                --outSAMmapqUnique {params.sam_mapq_unique} \
                --outSAMattributes All \
                --outSAMtype BAM Unsorted \
                --readFilesCommand zcat \
                --outReadsUnmapped Fastx \
                --outSAMattrRGline ID:{params.rg} SM:{params.rg} \
                --readFilesIn {input.fq1} > {log.f0}
            mv --force {params.prefix}Log.out {log.f1}
            mv --force {params.prefix}Log.final.out {output[2]}
            """

# ------------------------------------------------------------------------------
if config['generate_unmapped'] == True:
    rule gzip_unmapped_mate1:
        input:
            os.path.join(OD_UBAM, '{sample}_Unmapped.out.mate1'),
        output:
            os.path.join(OD, 'unmapped', '{sample}_1.fastq.gz')
        threads: 1           
        shell:
            """
            gzip -c {input} > {output}
            """

# ------------------------------------------------------------------------------
if config['generate_unmapped'] == True and config['library']['type'] == 'paired-end':
    rule gzip_unmapped_mate2:
        input:
            os.path.join(OD_UBAM, '{sample}_Unmapped.out.mate2'),
        output:
            os.path.join(OD, 'unmapped', '{sample}_2.fastq.gz')
        threads: 1
        shell:
            """
            gzip -c {input} > {output}
            """

# ------------------------------------------------------------------------------
"""
this sort rule is specific for STAR mapping (bulk rnaseq pipeline)
because of the bam file name: _Aligned.out.bam
"""
rule sortbam_star:
    input:
        os.path.join(OD_UBAM, '{sample}_Aligned.out.bam')
    output:
        (os.path.join(OD_BAM, '{sample}.bam'))
    log:
        os.path.join(OD_LOG, '{sample}.sortbam_star.log')
    threads: 1
    resources:
#        mem_mb = 120 * 1024
        mem_mb = 24000
    singularity:
        config['SAMTOOLS_IMAGE']
    shell:
        'rm -f {input}_sort.*.bam && \
          samtools sort -m 10G -@ {threads} -T {input}_sort -o {output} {input}'


# ------------------------------------------------------------------------------
# stats from STAR mapping
rule star_stats:
    input:
        pheno = {PHENODATA},
        files = expand(os.path.join(OD_LOG, '{sample}_Log.final.out'), sample=sample_ids)
    output:
         {MAPSTATS}
    threads: 1
    resources:
        mem_mb = 1000
    run:
        star_mapping_stats (output[0], input[0], OD_LOG)
