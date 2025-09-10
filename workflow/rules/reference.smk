# ------------------------------------------------------------------------------
# copy reference genome data from source to results directory
rule genome:
    input:
        gz = {GENOME_FASTA_GZ},
        fai = {GENOME_FAI},
        gzi = {GENOME_GZI},
        dict = {GENOME_DICT},
        sizes = {GENOME_SIZE},
    output:
        ugz = temp(os.path.join(OD_ANNO, 'genome.fa')),
        dict = temp(os.path.join(OD_ANNO, 'genome.fa.dict')),
        faiugz = temp(os.path.join(OD_ANNO, 'genome.fa.fai')),
        gz = os.path.join(OD_ANNO, 'genome.fa.gz'),
        fai = os.path.join(OD_ANNO, 'genome.fa.gz.fai'),
        gzi = os.path.join(OD_ANNO, 'genome.fa.gz.gzi'),
        sizes = os.path.join(OD_ANNO, 'genome.chrom.sizes'),
    params:
        custom_fasta = config.get("custom_fasta", None) if config.get("custom_fasta") else None  # Default to None if not set or empty
    threads: 1
    resources:
        mem_mb = 1000
    singularity:
        config['SAMTOOLS_IMAGE']      
    shell:
        """
        cp -L {input.gz} {output.gz}
        if [[ "{params.custom_fasta}" != "None" ]]; then
            cat {params.custom_fasta} | gzip -c >> {output.gz}
            gunzip -c {output.gz} > {output.ugz}
            samtools faidx {output.gz}    # gzi
            samtools faidx {output.ugz}   # fai
            samtools dict {output.ugz} > {output.dict}
            cut -f1,2 {output.ugz}.fai > {output.sizes}
        else
            cp -L {input.fai} {output.fai}
            cp -L {input.fai} {output.faiugz}
            cp -L {input.gzi} {output.gzi}
            cp -L {input.dict} {output.dict}
            cp -L {input.sizes} {output.sizes}
            gunzip -c {input.gz} > {output.ugz}
            samtools faidx {output.ugz} 
        fi
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
        dict = os.path.join(OD_ANNO, 'genome.fa.dict'),
    output:
        gtf = os.path.join(OD_ANNO, '{db}.gtf.gz'),
        ugtf = temp(os.path.join(OD_ANNO, '{db}.gtf')),
        len = os.path.join(OD_ANNO, '{db}.geneLength.gz'),
        flat = os.path.join(OD_ANNO, '{db}.refFlat.gz'),
        annot = os.path.join(OD_ANNO, '{db}.annot.gz'),
        loci = os.path.join(OD_ANNO, '{db}.loci.txt'),
        bed = os.path.join(OD_ANNO, '{db}.bed.gz'),
        ribo = os.path.join(OD_ANNO, '{db}.genome.rRNA_intervals'),
    params:
        custom_bed = config.get("custom_bed", None) if config.get("custom_bed") else None,
    threads: 1
    resources:
        mem_mb = 1000
    singularity:
        config['BEDOPS_IMAGE']
    shell:
        """
        cp -L {input.gtf} {output.gtf}
        cp -L {input.loci} {output.loci}
        cp -L {input.ribo} {output.ribo}
        cp -L {input.flat} {output.flat}
        gzip -c {input.len} > {output.len}
        gunzip -c {input.gtf} > {output.ugtf}
        gtf2bed < {output.ugtf} | gzip -c > {output.bed}
        grep -vw '^id' {input.annot} | awk 'BEGIN{{FS="\\t"; OFS="\\t"}} {{if (NF==4){{print $1, $2, $4}} else {{print $0}}}}' \
            | gzip -c > {output.annot}

        if [[ "{params.custom_bed}" != "None" ]]; then
            # BED file
            cat {params.custom_bed} | gzip -c >> {output.bed}

            # GTF file
            awk 'BEGIN {{FS="\\t"; OFS="\\t"}} $0 !~ /#/ \
                {{print $1, "BED_to_GTF", "gene", $2+1, $3, ".", $6, ".", "gene_id \\""$4"\\"; transcript_id \\""$4"\\"; gene_name \\""$4"\\";"}}' \
                {params.custom_bed} >> {output.ugtf}
            gzip -c {output.ugtf} > {output.gtf}

            # Gene lengths file
            awk 'BEGIN{{FS="\\t"; OFS="\\t"}} $0 !~ /#/ {{l=$3-$2+1; print $4, l, l, l, l}}' {params.custom_bed} \
                | gzip -c >> {output.len}

            # Flat file
            awk 'BEGIN {{FS="\\t"; OFS="\\t"}} $0 !~ /#/ {{print $4, $4, $1, $6, $2, $3, $2, $3, 1, $2, $3}}' {params.custom_bed} \
                | gzip -c >> {output.flat}

            # Annot file
            awk 'BEGIN{{FS="\\t"; OFS="\\t"}} $0 !~ /#/ {{print $4, $4, $4}}' {params.custom_bed} \
                | gzip -c >> {output.annot}
            
            # Loci file
            awk 'BEGIN{{FS="\\t"; OFS="\\t"}} $0 !~ /#/ {{print $1, $2, $3, $6, $4, $4, $4}}' {params.custom_bed} >> {output.loci}

            # Ribo file (modify only header)
            grep '^@SQ' {input.dict} | cut -f1-3 > {output.ribo}
            grep -v '^@SQ' {input.ribo} >> {output.ribo}
        fi
        """

# ------------------------------------------------------------------------------
# Rule to generate STAR index from the copied reference genome
rule star_index:
    input:
        genome_fasta = rules.genome.output.ugz,
        chrom_sizes = rules.genome.output.sizes,
    output:
        genome_dir = protected(directory(os.path.join(OD_ANNO, "STAR_Index"))),
    params:
        tmp_dir = os.path.join(OD_ANNO, "__STARtmp"),
    threads: config['star_index_threads']
    resources:
        mem_mb = config['star_index_mem_mb']
    singularity:
        config['STAR_IMAGE']
    log:
        log_file = os.path.join(OD_LOG, "star_index.log"),
    shell:
        """
        # Calculate genomeSAindexNbases dynamically
        GENOME_SA_INDEX_NBASES=$(awk -F '\\t' 'NR==1{{TOT=1}} NR>1{{TOT+=$2}} END {{VAL=int(((log(TOT)/log(2))/2)-1); if (VAL < 14) {{print VAL}} else {{print 14}}}}' {input.chrom_sizes})

        # Build the STAR index
        STAR --runMode genomeGenerate \
             --runThreadN {threads} \
             --outTmpDir {params.tmp_dir} \
             --genomeDir {output.genome_dir} \
             --genomeSAindexNbases $GENOME_SA_INDEX_NBASES \
             --genomeFastaFiles {input.genome_fasta} \
             2>&1 | tee {log.log_file}
       """