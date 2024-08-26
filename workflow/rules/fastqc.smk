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
if config['library']['type'] == 'paired-end':
    rule validate_fastq:
        input:
            os.path.join(OD_FASTQ, '{name}_1.fastq.gz'),
            os.path.join(OD_FASTQ, '{name}_2.fastq.gz'),
            os.path.join(OD_FASTQ, '{name}_1.fastq.gz.validated'),
            os.path.join(OD_FASTQ, '{name}_2.fastq.gz.validated')
        output:
            temp(os.path.join(OD_FASTQ, '{name}_1.fastq.gz.done')),
            temp(os.path.join(OD_FASTQ, '{name}_2.fastq.gz.done'))
        threads: 1
        resources:
            mem_mb=1000
        run:
            # check if first from both mate files have same name
            import gzip
            import re

            with gzip.open(input[0],'r') as fin: 
                first_line1 = fin.readline().strip()
            with gzip.open(input[1],'r') as fin: 
                first_line2 = fin.readline().strip()
            r1 = re.split(' ', first_line1.decode('utf-8'))[0]
            r2 = re.split(' ', first_line2.decode('utf-8'))[0]
            if r1 != r2:
                print('ERROR: First read names are different!', file=sys.stderr)
                print('       read 1:'+r1, file=sys.stderr)
                print('       read 2:'+r2, file=sys.stderr)

            n0 = count_gzip_lines(input[0])
            n1 = count_gzip_lines(input[1])
            if n0 == n1:
                if n0 % 4 == 0:
                    Path(output[0]).touch()
                    Path(output[1]).touch()
                else:
                    print('ERROR: Number of lines in fastq file is not multiple of four: ', n0, file=sys.stderr)
            else:
                print('ERROR: Number of reads in mate files are different:', file=sys.stderr)
                print('File '+input[0]+' contains '+str(n0)+' lines', file=sys.stderr)
                print('File '+input[1]+' contains '+str(n1)+' lines', file=sys.stderr)
else:
    rule validate_fastq:
        input:
            os.path.join(OD_FASTQ, '{name}.fastq.gz'), 
            os.path.join(OD_FASTQ, '{name}.fastq.gz.validated')
        output:
            temp(os.path.join(OD_FASTQ, '{name}.fastq.gz.done'))
        threads: 1
        resources:
            mem_mb=1000
        run:
            n0 = count_gzip_lines(input[0])
            if n0 % 4 == 0:
                Path(output[0]).touch()
            else:
                print('ERROR: Number of lines in fastq file is not multiple of four:', n0, file=sys.stderr)


# ---------------------------------------------------------------
if config['cutadapt']['run'] == True:
    rule cutadapt:
        input:
            os.path.join(OD_FASTQ, '{name}_1.fastq.gz'),
            os.path.join(OD_FASTQ, '{name}_2.fastq.gz')
        output:
            fq1 = temp(os.path.join(OD_CUTADAPT, '{name}_1.fastq.gz')),
            fq2 = temp(os.path.join(OD_CUTADAPT, '{name}_2.fastq.gz')),
            report = temp(os.path.join(OD_CUTADAPT, '{name}.report.txt')),
        log:
            os.path.join(OD_LOG, '{name}.cutadapt.log')
        threads: 12
        resources:
            mem_mb=10000
        params:
            config['cutadapt']['parameters']
        singularity:
            config['CUTADAPT_IMAGE']
        shell:
            """
            cutadapt --cores {threads} {params} -o {output[0]} -p {output[1]} {input[0]} {input[1]} > {output.report}
            """
else:
    if config['library']['type'] == 'paired-end':
        rule no_cutadapt:
            input:
                os.path.join(OD_FASTQ, '{name}_1.fastq.gz'),
                os.path.join(OD_FASTQ, '{name}_2.fastq.gz')
            output:
                report = temp(os.path.join(OD_CUTADAPT, '{name}.report.txt'))
            shell:
                """
                touch {output}
                """
    else:
        rule no_cutadapt:
            input:
                os.path.join(OD_FASTQ, '{name}.fastq.gz')
            output:
                report = temp(os.path.join(OD_CUTADAPT, '{name}.report.txt'))
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
