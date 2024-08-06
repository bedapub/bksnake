# files = lambda wildcards: expand(os.path.join(OD_FC, '{sample}.{db}.cnt.gz'), sample=sample_ids, db=wildcards.db),

# -----------------------------------------------------------------------------------------
if config['pipeline'] == 'bksnake':
    rule multiqc: # Rule for RefSeq/Ensembl annotations
        input:
            cnts = lambda wildcards: expand(os.path.join(OD_FC, '{sample}.{db}.cnt.gz'), sample=sample_ids, db=wildcards.db),
            summaries = lambda wildcards: expand(os.path.join(OD_FC, '{sample}.{db}.cnt.summary'), sample=sample_ids, db=wildcards.db),
            metrics = lambda wildcards: expand(os.path.join(OD_METRICS, '{sample}.{db}.RNAmetrics.txt'), sample=sample_ids, db=wildcards.db),
            rseqc = expand(os.path.join(OD_METRICS, '{sample}.strandedness.txt'), sample=sample_ids),
            html = expand(os.path.join(OD_FASTQC, '{name}_fastqc.html'), name=fastq_names_noext),
            fastqc = expand(os.path.join(OD_FASTQC, '{name}_fastqc'), name=fastq_names_noext),
            flagstat = expand(os.path.join(OD_STATS, '{sample}.bam.flagstat'), sample=sample_ids),
            stats = expand(os.path.join(OD_STATS, '{sample}.bam.stats'), sample=sample_ids),
            stats2 = expand(os.path.join(OD_STATS, '{sample}.bam.stats2'), sample=sample_ids),
            final = expand(os.path.join(OD_LOG, '{sample}_Log.final.out'), sample=sample_ids),
            cutadapt = expand(os.path.join(OD_CUTADAPT, '{sample}.report.txt'), sample=sample_ids),
        output:
            html = report(os.path.join(OD, '{db}_multiqc_report.html'), 
                htmlindex=os.path.join(OD, '{db}_multiqc_report.html'), 
                category='Quality Control', subcategory='{db}', caption='../report/fig3_multiqc.rst',
                labels={
                    'annotation': '{db}',
                    'figure': 'MultiQC Report'
                    }),
            metrics = os.path.join(OD_MULTIQC, '{db}/multiqc_picard_RnaSeqMetrics.txt'),
            counts = os.path.join(OD_MULTIQC, '{db}/multiqc_featureCounts.txt'),
            star = os.path.join(OD_MULTIQC, '{db}/multiqc_star.txt'),
            fastqc = os.path.join(OD_MULTIQC, '{db}/multiqc_fastqc.txt'),
        log:
            os.path.join(OD_LOG, 'multiqc_{db}.log')
        threads: 1
        resources:
            mem_mb=20000
        params:
            html = '{db}_multiqc_report.html',
            outdir = os.path.join(OD_MULTIQC, '{db}'),
            title = 'Report using {db} annotations',
            comment = 'MultiQC report utilized only {db} gene annotations.',
        singularity:
            config['MULTIQC_IMAGE']
        shell:
            """
            #unset LD_PRELOAD
            LC_ALL='en_US.utf-8'
            multiqc --no-megaqc-upload --verbose --force \
                --title '{params.title}' \
                --comment '{params.comment}' \
                --outdir {params.outdir} \
                --filename {params.html} \
                {OD}/cutadapt \
                {OD}/fastqc \
                {OD}/log/*_Log.final.out \
                {OD}/metrics \
                {OD}/fc \
                {input.summaries} 2> {log} \
                && mkdir -p {params.outdir} \
                && mv -f {params.outdir}/{wildcards.db}_multiqc_report.html {output.html} \
                && mv -f {params.outdir}/{wildcards.db}_multiqc_report_data/* {params.outdir}/ \
                && rmdir {params.outdir}/{wildcards.db}_multiqc_report_data
            """

elif config['pipeline'] == 'vcsnake':
    rule multiqc:
        input:
            fastqc_html = expand(os.path.join(OD_FASTQC, '{name}_fastqc.html'), name=fastq_names_noext),
            vep_html = expand(os.path.join(OD_VCF, '{sample}.vep.all_chroms_summary.html'), sample=sample_ids)
        output:
            os.path.join(OD, 'multiqc_report.html')
        log:
            os.path.join(OD_LOG, 'multiqc.log')
        threads: 1
        resources:
            mem_mb=20000
        params:
            report = 'multiqc_report.html',
            ignores = '-x '+OD_BW+' -x unmapped -x bam '\
                      '-x '+OD_ANNO+' -x jbrowse -x fastq'
        singularity:
            config['MULTIQC_IMAGE']
        shell:
            """
            #unset LD_PRELOAD
            LC_ALL='en_US.utf-8'
            multiqc --no-megaqc-upload --verbose {params.ignores} --force -o {OD} {OD} 2> {log}
            """
else:
    sys.exit('Variable \'pipeline\' not defined in input config yaml file.')  
