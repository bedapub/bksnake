# copy reference genome data from source to results directory
rule refdata:
    input:
        gtf_ref = {GTF_REF},
        gtf_ens = {GTF_ENS},
        len_ref = {LEN_REF},
        len_ens = {LEN_ENS},
        annot_ref = {ANNOT_REF},
        annot_ens = {ANNOT_ENS},
        loci_ref = {LOCI_REF},
        loci_ens = {LOCI_ENS},
        gz = {GENOME_FASTA_GZ},
        fai = {GENOME_FAI},
        gzi = {GENOME_GZI},
        flat_ref = {REFFLAT_REF},
        flat_ens = {REFFLAT_ENS},
        ribo = {RIBO_INTERVALS}
    output:
        gtf_ref = os.path.join(OD_ANNO,'refseq.gtf.gz'),
        gtf_ens = os.path.join(OD_ANNO,'ensembl.gtf.gz'),
        len_ref = os.path.join(OD_ANNO,'refseq.geneLength.gz'),
        len_ens = os.path.join(OD_ANNO,'ensembl.geneLength.gz'),
        annot_ref = os.path.join(OD_ANNO,'refseq.annot.gz'),
        annot_ens = os.path.join(OD_ANNO,'ensembl.annot.gz'),
        loci_ref = os.path.join(OD_ANNO,'refseq.loci.txt'),
        loci_ens = os.path.join(OD_ANNO,'ensembl.loci.txt'),
        gz = os.path.join(OD_ANNO,'genome.fa.gz'),
        fai = os.path.join(OD_ANNO,'genome.fa.gz.fai'),
        gzi = os.path.join(OD_ANNO,'genome.fa.gz.gzi'),
        flat_ref = os.path.join(OD_ANNO,'refseq.refFlat.gz'),
        flat_ens = os.path.join(OD_ANNO,'ensembl.refFlat.gz'),
        ribo = os.path.join(OD_ANNO,'genome.rRNA_intervals')
    threads: 1
    resources:
        mem_mb = 1000
    shell:
        """
        gzip -c {input.gtf_ref} > {output.gtf_ref}
        gzip -c {input.gtf_ens} > {output.gtf_ens}
        gzip -c {input.len_ref} > {output.len_ref}
        gzip -c {input.len_ens} > {output.len_ens}
        gzip -c {input.annot_ref} > {output.annot_ref}
        gzip -c {input.annot_ens} > {output.annot_ens}
        cp -pr {input.loci_ref} {output.loci_ref}
        cp -pr {input.loci_ens} {output.loci_ens}
        cp -pr {input.gz} {output.gz}
        cp -pr {input.fai} {output.fai}
        cp -pr {input.gzi} {output.gzi}
        gzip -c {input.flat_ref} > {output.flat_ref}
        gzip -c {input.flat_ens} > {output.flat_ens}
        cp -pr {input.ribo} {output.ribo}
        """
