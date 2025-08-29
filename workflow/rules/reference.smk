# ------------------------------------------------------------------------------
# copy reference genome data from source to results directory
rule genome:
    input:
        gz = {GENOME_FASTA_GZ},
        fai = {GENOME_FAI},
        gzi = {GENOME_GZI},
        dict = {GENOME_DICT},
    output:
        ugz = temp(os.path.join(OD_ANNO, 'genome.fa')),
        dict = temp(os.path.join(OD_ANNO, 'genome.fa.dict')),
        faiugz = temp(os.path.join(OD_ANNO, 'genome.fa.fai')),
        gz = os.path.join(OD_ANNO, 'genome.fa.gz'),
        fai = os.path.join(OD_ANNO, 'genome.fa.gz.fai'),
        gzi = os.path.join(OD_ANNO, 'genome.fa.gz.gzi'),
    threads: 1
    resources:
        mem_mb = 1000
    singularity:
        config['SAMTOOLS_IMAGE']      
    shell:
        """
        cp -L {input.gz} {output.gz}
        cp -L {input.fai} {output.fai}
        cp -L {input.fai} {output.faiugz}
        cp -L {input.gzi} {output.gzi}
        cp -L {input.dict} {output.dict}
        gunzip -c {input.gz} > {output.ugz}
        samtools faidx {output.ugz} 
        """

# ------------------------------------------------------------------------------
# copy reference genome annotations from source to results directory
rule annotations:
    input:
        gtf = os.path.join(GENOME_DIR, '{db}/annotation.gtf.gz'),
        len = os.path.join(GENOME_DIR, '{db}/annotation.geneLength'),
        flat = os.path.join(GENOME_DIR, '{db}/annotation.refFlat.gz'),
        annot = os.path.join(GENOME_DIR, '{db}/genes.tsv'),
        loci = os.path.join(GENOME_DIR, '{db}/genes.loci.txt'),
        ribo = os.path.join(GENOME_DIR, '{db}/annotation.rRNA.interval_list'),
    output:
        gtf = os.path.join(OD_ANNO, '{db}.gtf.gz'),
        ugtf = temp(os.path.join(OD_ANNO, '{db}.gtf')),
        len = os.path.join(OD_ANNO, '{db}.geneLength.gz'),
        flat = os.path.join(OD_ANNO, '{db}.refFlat.gz'),
        annot = os.path.join(OD_ANNO, '{db}.annot.gz'),
        loci = os.path.join(OD_ANNO, '{db}.loci.txt'),
        bed = os.path.join(OD_ANNO, '{db}.bed.gz'),
        ribo = os.path.join(OD_ANNO, '{db}.genome.rRNA_intervals'),
    threads: 1
    resources:
        mem_mb = 1000
    singularity:
        config['BEDOPS_IMAGE']
    shell:
        """
        cp -L {input.gtf} {output.gtf}
        gunzip -c {input.gtf} > {output.ugtf}
        gzip -c {input.len} > {output.len}
        cp -L {input.flat} {output.flat}
        grep -vw '^id' {input.annot} | awk 'BEGIN{{FS=\"\\t\"}}{{if (NF==4){{printf(\"%s\\t%s\\t%s\\n\", $1,$2,$4)}}else{{printf(\"%s\\n\",$0)}}}}' | gzip -c > {output.annot}
        cp -L {input.loci} {output.loci}
        gtf2bed < {output.ugtf} | gzip -c > {output.bed}
        cp -L {input.ribo} {output.ribo}
        """
