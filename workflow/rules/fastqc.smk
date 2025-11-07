"""
QC analysis for fastq files
"""

"""------------------------------------------------------------------------------
Define helper functions
For counting number of lines in (fastq) file
"""
def count_gzip_lines(filename):
    """Counts the number of lines in an gzipped input file.

    Parameters
    ----------
    filename : str
        Input file name.
        
    Returns
    -------
    i : int
        Number of lines in input file.
    """
    import gzip

    if not os.path.isfile(filename):
        raise Exception('count_gzip_lines: file is not present:'+filename)
    elif os.path.getsize(filename) == 0:
        raise Exception('count_gzip_lines: file is empty, i.e. getsize is zero:'+filename)
    else:
        with gzip.open(filename, 'rb') as f:
            for i, l in enumerate(f):
                pass
        return i+1


# ---------------------------------------------------------------
# determine md5 check sum for each fastq file
rule md5sum:
    input:
        os.path.join(OD_FASTQ, '{name}')
    output:
        temp(os.path.join(OD, 'md5sum', '{name}.md5sum'))
    threads: 1
    resources:
        mem_mb=1000
    shell:
        """
        md5sum {input} > {output}
        """

# ---------------------------------------------------------------
# compare the md5 files
rule compare_md5sum:
    input:
        files = expand(os.path.join(OD, 'md5sum', '{name}.fastq.gz.md5sum'), name=fastq_names_noext)
    output:
        md5sum = os.path.join(OD_QC, 'md5sum.txt')
    params:
        dir = os.path.join(OD, 'md5sum')
    shell:
        """
        cat {input.files} > {output}
        workflow/scripts/md5sum.sh {params.dir} >> {output}
        """

# ---------------------------------------------------------------
# check integrity of raw input gzipped fastq files
rule validate_gzip:
    input:
        os.path.join(OD_FASTQ, '{name}')
    output:
        temp(os.path.join(OD_FASTQ, '{name}.validated'))
    threads: 1
    resources:
        mem_mb=1000
    shell:
        """
        gzip -t {input} && touch {output}
        """

# ---------------------------------------------------------------
rule validate_fastq_paired:
    input:
        os.path.join(OD_FASTQ, '{name}_1.fastq.gz'),
        os.path.join(OD_FASTQ, '{name}_2.fastq.gz'),
        os.path.join(OD_FASTQ, '{name}_1.fastq.gz.validated'),
        os.path.join(OD_FASTQ, '{name}_2.fastq.gz.validated'),
    output:
        temp(os.path.join(OD_FASTQ, '{name}_1.fastq.gz.done')),
        temp(os.path.join(OD_FASTQ, '{name}_2.fastq.gz.done')),
    wildcard_constraints:
        name = "|".join(SAMPLES_PE) if SAMPLES_PE else "DUMMY_NEVER_MATCH"
    threads: 1
    resources:
        mem_mb=1000
    shell:
        """
        python3 workflow/scripts/validate_fastq.py paired {input[0]} {input[1]} {output[0]} {output[1]}
        """

# ---------------------------------------------------------------
rule validate_fastq_single:
    input:
        os.path.join(OD_FASTQ, '{name}.fastq.gz'),
        os.path.join(OD_FASTQ, '{name}.fastq.gz.validated'),
    output:
        temp(os.path.join(OD_FASTQ, '{name}.fastq.gz.done')),
    wildcard_constraints:
        name = "|".join(SAMPLES_SE) if SAMPLES_SE else "DUMMY_NEVER_MATCH"
    threads: 1
    resources:
        mem_mb=1000
    shell:
        """
        python3 workflow/scripts/validate_fastq.py single {input[0]} {output[0]}
        """

# ---------------------------------------------------------------
if config['cutadapt']['run'] == True:
    rule cutadapt_paired:
        input:
            os.path.join(OD_FASTQ, '{name}_1.fastq.gz'),
            os.path.join(OD_FASTQ, '{name}_2.fastq.gz'),
        output:
            fq1 = temp(os.path.join(OD_CUTADAPT, '{name}_1.fastq.gz')),
            fq2 = temp(os.path.join(OD_CUTADAPT, '{name}_2.fastq.gz')),
            report = temp(os.path.join(OD_CUTADAPT, '{name}.report.txt')),
        wildcard_constraints:
            name = "|".join(SAMPLES_PE) if SAMPLES_PE else "DUMMY_NEVER_MATCH"
        log:
            os.path.join(OD_LOG, '{name}.cutadapt_paired.log')
        threads: 12
        resources:
            mem_mb=10000
        params:
            config['cutadapt']['parameters_paired']
        singularity:
            config['CUTADAPT_IMAGE']
        shell:
            """
            cutadapt --cores {threads} {params} -o {output[0]} -p {output[1]} {input[0]} {input[1]} > {output.report}
            """

    rule cutadapt_single:
        input:
            fq = os.path.join(OD_FASTQ, '{name}.fastq.gz'),
        output:
            fq = temp(os.path.join(OD_CUTADAPT, '{name}.fastq.gz')),
            report = temp(os.path.join(OD_CUTADAPT, '{name}.report.txt')),
        wildcard_constraints:
            name = "|".join(SAMPLES_SE) if SAMPLES_SE else "DUMMY_NEVER_MATCH"
        log:
            os.path.join(OD_LOG, '{name}.cutadapt_single.log')
        threads: 12
        resources:
            mem_mb=10000
        params:
            config['cutadapt']['parameters_single']
        singularity:
            config['CUTADAPT_IMAGE']
        shell:
            """
            cutadapt --cores {threads} {params} -o {output.fq} {input.fq} > {output.report}
            """
else:
    rule no_cutadapt_paired:
        input:
            os.path.join(OD_FASTQ, '{name}_1.fastq.gz'),
            os.path.join(OD_FASTQ, '{name}_2.fastq.gz')
        output:
            report = temp(os.path.join(OD_CUTADAPT, '{name}.report.txt'))
        wildcard_constraints:
            name = "|".join(SAMPLES_PE) if SAMPLES_PE else "DUMMY_NEVER_MATCH"
        shell:
            """
            touch {output}
            """
    rule no_cutadapt_single:
        input:
            os.path.join(OD_FASTQ, '{name}.fastq.gz')
        output:
            report = temp(os.path.join(OD_CUTADAPT, '{name}.report.txt'))
        wildcard_constraints:
            name = "|".join(SAMPLES_SE) if SAMPLES_SE else "DUMMY_NEVER_MATCH"
        shell:
            """
            touch {output}
            """

# ---------------------------------------------------------------
rule fastqc:
    input:      
        os.path.join(OD_CUTADAPT, '{name}.fastq.gz') if config['cutadapt']['run'] else os.path.join(OD_FASTQ, '{name}.fastq.gz')
    output:
        html = os.path.join(OD_FASTQC, '{name}_fastqc.html'),
        dir = temp(directory(os.path.join(OD_FASTQC, '{name}_fastqc'))),
        zip = temp(os.path.join(OD_FASTQC, '{name}_fastqc.zip'))
    log:
        os.path.join(OD_LOG, '{name}.fastqc.log')
    threads: 4
    resources:
        mem_mb=10000
    params:
       outdir = OD_FASTQC
    singularity:
        config['FASTQC_IMAGE']
    shell:
        """
        fastqc --extract --threads {threads} \
            --outdir {params.outdir} {input[0]} 2> {log}
        """
