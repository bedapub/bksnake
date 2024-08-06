# ------------------------------------------------------------------------------
# copy reference genome data from source to results directory
rule genome:
    input:
        gz = {GENOME_FASTA_GZ},
        fai = {GENOME_FAI},
        gzi = {GENOME_GZI},
        ribo = {RIBO_INTERVALS}
    output:
        gz = os.path.join(OD_ANNO, 'genome.fa.gz'),
        fai = os.path.join(OD_ANNO, 'genome.fa.gz.fai'),
        gzi = os.path.join(OD_ANNO, 'genome.fa.gz.gzi'),
        ribo = os.path.join(OD_ANNO, 'genome.rRNA_intervals')
    threads: 1
    resources:
        mem_mb = 1000
    shell:
        """
        cp -Lpr {input.gz} {output.gz}
        cp -Lpr {input.fai} {output.fai}
        cp -Lpr {input.gzi} {output.gzi}
        cp -Lpr {input.ribo} {output.ribo}
        """

# ------------------------------------------------------------------------------
# copy reference genome annotations from source to results directory
rule annotations:
    input:
        gtf = os.path.join(GENOME_DIR, 'gtf/{db}/{db}.exons.gtf'),
        len = os.path.join(GENOME_DIR, 'gtf/{db}/{db}.exons.geneLength'),
        flat = os.path.join(GENOME_DIR, 'gtf/{db}/{db}.refFlat'),
        annot = os.path.join(GENOME_DIR, 'gtf/{db}/{db}.exons.annot'),
        loci = os.path.join(GENOME_DIR, 'gtf/{db}/{db}.loci.txt'),
    output:
        gtf = os.path.join(OD_ANNO, '{db}.gtf.gz'),
        len = os.path.join(OD_ANNO, '{db}.geneLength.gz'),
        flat = os.path.join(OD_ANNO, '{db}.refFlat.gz'),
        annot = os.path.join(OD_ANNO, '{db}.annot.gz'),
        loci = os.path.join(OD_ANNO, '{db}.loci.txt'),
        bed = os.path.join(OD_ANNO, '{db}.bed'),
    threads: 1
    resources:
        mem_mb = 1000
    singularity:
        config['BEDOPS_IMAGE']
    shell:
        """
        gzip -c {input.gtf} > {output.gtf}
        gzip -c {input.len} > {output.len}
        gzip -c {input.flat} > {output.flat}
        gzip -c {input.annot} > {output.annot}
        cp -Lpr {input.loci} {output.loci}        
        gtf2bed < {input.gtf} > {output.bed}
        """
