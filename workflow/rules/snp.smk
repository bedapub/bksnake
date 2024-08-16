"""
Rule for cross-checking sample matching using SNP fingerprints (snpMatch) 

and summarizing results from multiple checks (concatMetrics)


Source: https://github.com/naumanjaved/fingerprint_maps

Citation: 
1. Javed, N., Farjoun, Y., Fennell, T.J. et al. Detecting sample swaps in diverse NGS data types using linkage disequilibrium. Nat Commun 11, 3697 (2020). DOI: https://doi.org/10.1038/s41467-020-17453-5
2. Picard Tool from the Broad Institute

How to make the haplotype map:

  wget -O hg38_chr.map  https://raw.githubusercontent.com/naumanjaved/fingerprint_maps/master/map_files/hg38_chr.map
  cat <genome directory>/hg38/fasta/genome.dict > resources/haplotype.map
  grep -v '^@' hg38_chr.map >> resources/haplotype.map

"""

MAP = os.path.join(OD, 'haplotype.map')

"""
Create the map file
"""
rule haplotype_map:
    output:
        tmp = temp(os.path.join(OD, 'hg38_chr.map')),
        map = temp(MAP)
    params:
        uri = 'https://raw.githubusercontent.com/naumanjaved/fingerprint_maps/master/map_files/hg38_chr.map'
    shell:
        """
        wget -q -N -O {output.tmp} {params.uri} \
            && cat {GENOME_DICT} > {output.map} \
            && grep -v '^@' {output.tmp} >> {output.map}
        """


"""
This rule is not required because the RG is added by STAR
"""
rule addreplacerg:
    input:        
        bam = os.path.join(OD_BAM, '{sample}.bam'),
        bai = os.path.join(OD_BAM, '{sample}.bam.bai')
    output:
        bam = temp(os.path.join(OD_VCF, '{sample}.rg.bam')),
        bai = temp(os.path.join(OD_VCF, '{sample}.rg.bam.bai'))
    log:
        os.path.join(OD_LOG, '{sample}.addreplacerg.log')
    params:
        id = '{sample}'
    threads: 8
    resources:
        mem_mb = 20000
    singularity:
        config['SAMTOOLS_IMAGE']
    shell:
        """
        samtools addreplacerg --threads {threads} -o {output.bam} -r '@RG\\tID:{params.id}\\tSM:{params.id}\\tLB:{params.id}\\tPL:Illumina' {input.bam} \
        && samtools index --threads {threads} {output.bam}
        """

"""
Determine finger prints by picard
"""
rule fingerprint:
    input:
        map = MAP,
        ref = rules.genome.output.ugz,
        fai = rules.genome.output.faiugz,
        dict = rules.genome.output.dict,
        bam = os.path.join(OD_BAM, '{sample}.bam'),
        bai = os.path.join(OD_BAM, '{sample}.bam.bai')
    output:
        vcf = os.path.join(OD_VCF, '{sample}.marked_duplicates_fingerprint.vcf') if config['keep_vcf_files'] else temp(os.path.join(OD_VCF, '{sample}.marked_duplicates_fingerprint.vcf')),
        idx = os.path.join(OD_VCF, '{sample}.marked_duplicates_fingerprint.vcf.idx') if config['keep_vcf_files'] else temp(os.path.join(OD_VCF, '{sample}.marked_duplicates_fingerprint.vcf.idx')),
    log:
        os.path.join(OD_LOG, '{sample}.fingerprint.log')
    params:
        map = MAP
    threads: 1
    resources:
        mem_mb = 20000
    singularity:
        config['PICARD_IMAGE']
    shell:
        """
        export _JAVA_OPTIONS="-Xmx10g" && \
        /usr/local/bin/picard ExtractFingerprint \
            -INPUT {input.bam} \
            -REFERENCE_SEQUENCE {input.ref} \
            -HAPLOTYPE_MAP {input.map} \
            -OUTPUT {output.vcf} 2> {log}
        """


"""
Crosscheck fingerprints
"""
def convert_to_param(input_files_list):
    files = ' '.join(input_files_list).split()
    files = ['--INPUT ' + file for file in files]
    return ' '.join(files)

rule crosscheck:
    input:
        map = MAP,
        ref = rules.genome.output.ugz,
        fai = rules.genome.output.faiugz,
        dict = rules.genome.output.dict,
        vcfs = expand(os.path.join(OD_VCF, '{sample}.marked_duplicates_fingerprint.vcf'), sample=sample_ids),
        idxs = expand(os.path.join(OD_VCF, '{sample}.marked_duplicates_fingerprint.vcf.idx'), sample=sample_ids)
    output:
        file = os.path.join(OD_METRICS, 'crosscheck_metrics')
    log:
        os.path.join(OD_LOG, 'crosscheck.log')
    params:
        type = 'SAMPLE',
        input = convert_to_param(expand(os.path.join(OD_VCF, '{sample}.marked_duplicates_fingerprint.vcf'), sample=sample_ids))
    threads: 4
    resources:
        mem_mb = 20000
    singularity:
        config['PICARD_IMAGE']
    shell:
        """
        export _JAVA_OPTIONS="-Xmx10g" \
        && /usr/local/bin/picard CrosscheckFingerprints \
           --EXIT_CODE_WHEN_MISMATCH 0 \
           --EXIT_CODE_WHEN_NO_VALID_CHECKS 0 \
           --HAPLOTYPE_MAP {input.map} \
           --REFERENCE_SEQUENCE {input.ref} \
           --CROSSCHECK_BY {params.type} \
           --NUM_THREADS {threads} \
           --OUTPUT {output.file} \
           {params.input} 2> {log}
        """

