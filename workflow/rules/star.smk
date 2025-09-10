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
# Try calculate memory for star per thread: we want in total 120Gb
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html
def get_mem_mb(wildcards, threads):
    return 2000 * threads # 20G for each threads
    #return 120 * 1024 / threads # 120G / threads -->FAILS


def star_input_done_files(wildcards):
    files = []
    if config['library']['type'] == 'paired-end':
        files.append(os.path.join(OD_FASTQ, f"{wildcards.sample}_1.fastq.gz.done"))
        files.append(os.path.join(OD_FASTQ, f"{wildcards.sample}_2.fastq.gz.done"))
    else:
        files.append(os.path.join(OD_FASTQ, f"{wildcards.sample}.fastq.gz.done"))
    return files


def star_input_fastq_files(wildcards):
    files = []
    if config['library']['type'] == 'paired-end':
        if config['cutadapt']['run']:
            files.append(os.path.join(OD_CUTADAPT, f"{wildcards.sample}_1.fastq.gz"))
            files.append(os.path.join(OD_CUTADAPT, f"{wildcards.sample}_2.fastq.gz"))
        else:
            files.append(os.path.join(OD_FASTQ, f"{wildcards.sample}_1.fastq.gz"))
            files.append(os.path.join(OD_FASTQ, f"{wildcards.sample}_2.fastq.gz"))
    else:
        if config['cutadapt']['run']:
            files.append(os.path.join(OD_CUTADAPT, f"{wildcards.sample}.fastq.gz"))
        else:
            files.append(os.path.join(OD_FASTQ, f"{wildcards.sample}.fastq.gz"))
    return files


# ------------------------------------------------------------------------------
rule star:
    input:
        fq = star_input_fastq_files,
        done = star_input_done_files,
        gtfs = expand(rules.annotations.output.ugtf, db=DBS),
        cutadapt = os.path.join(OD_CUTADAPT, '{sample}.report.txt'),
        genome_dir = rules.star_index.output.genome_dir if config.get('generate_star_index', False) else STAR_DIR,
    output:
        tmp1 = temp(directory(os.path.join(OD_UBAM, '{sample}__STARtmp'))),
        tmp2 = temp(directory(os.path.join(OD_UBAM, '{sample}__STARgenome'))),
        out = temp(os.path.join(OD_LOG, '{sample}_Log.final.out')),
        bam = temp(os.path.join(OD_UBAM, '{sample}_Aligned.out.bam')),
        log = temp(os.path.join(OD_UBAM, '{sample}_Log.progress.out')),
        tab = temp(os.path.join(OD_UBAM, '{sample}_SJ.out.tab')),
        mate1 = temp(os.path.join(OD_UBAM, '{sample}_Unmapped.out.mate1')),
        mate2 = temp(os.path.join(OD_UBAM, '{sample}_Unmapped.out.mate2')) if config['library']['type'] == 'paired-end' else [],            
    log:
        f0 = os.path.join(OD_LOG, '{sample}.star.log'),
        f1 = os.path.join(OD_LOG, '{sample}_Log.out')
    retries: 3
    threads: config['star_threads']
    resources:
        mem_mb = get_mem_mb
    params:
        star_params = config['star']['params'],
        prefix = join(OD_UBAM, '{sample}_'),
        rg = '{sample}'
    singularity:
        config['STAR_IMAGE']
    shell:
        """
        STAR --runThreadN {threads} \
            {params.star_params} \
            --genomeDir {input.genome_dir} \
            --outFileNamePrefix {params.prefix} \
            --outTmpDir {output.tmp1} \
            --sjdbGTFfile {GTF_FOR_STAR_MAPPING} \
            --outSAMtype BAM Unsorted \
            --readFilesCommand zcat \
            --outReadsUnmapped Fastx \
            --outSAMattrRGline ID:{params.rg} SM:{params.rg} \
            --readFilesIn {input.fq} > {log.f0}
        mv --force {params.prefix}Log.out {log.f1}
        mv --force {params.prefix}Log.final.out {output.out}
        """

# ------------------------------------------------------------------------------
if config['generate_unmapped'] == True:
    rule gzip_unmapped_mate1:
        input:
            rules.star.output.mate1,
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
            rules.star.output.mate2,
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
        rules.star.output.bam,
    output:
        os.path.join(OD_BAM, '{sample}.bam') if config['keep_bam_files'] else temp(os.path.join(OD_BAM, '{sample}.bam')),
    log:
        os.path.join(OD_LOG, '{sample}.sortbam_star.log')
    threads: 1
    resources:
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
        files = expand(rules.star.output.out, sample=sample_ids),
    output:
         {MAPSTATS}
    threads: 1
    resources:
        mem_mb = 1000
    run:
        star_mapping_stats (output[0], input[0], OD_LOG)
