
# -----------------------------------------------------------------------------------------
if config['pipeline'] == 'bksnake':
    rule multiqc_ref: # Rule for RefSeq annotations
        input:
            expand(os.path.join(OD_METRICS, '{sample}.refseq.RNAmetrics.txt'), sample=sample_ids),
            expand(os.path.join(OD_FASTQC, '{name}_fastqc.html'), name=fastq_names_noext),
            expand(os.path.join(OD_FASTQC, '{name}_fastqc'), name=fastq_names_noext),
            expand(os.path.join(OD_STATS, '{sample}.bam.flagstat'), sample=sample_ids),
            expand(os.path.join(OD_STATS, '{sample}.bam.stats'), sample=sample_ids),
            expand(os.path.join(OD_STATS, '{sample}.bam.stats2'), sample=sample_ids),
            expand(os.path.join(OD_LOG, '{sample}_Log.final.out'), sample=sample_ids),
            expand(os.path.join(OD_FC, '{sample}.refseq.cnt.summary'), sample=sample_ids),
            expand(os.path.join(OD_FC, '{sample}.refseq.cnt.gz'), sample=sample_ids),
        output:
            os.path.join(OD, 'multiqc_report.html'),
            os.path.join(OD_MULTIQC, 'multiqc_picard_RnaSeqMetrics.txt'),
            os.path.join(OD_MULTIQC, 'multiqc_featureCounts.txt'),
            os.path.join(OD_MULTIQC, 'multiqc_star.txt'),
            os.path.join(OD_MULTIQC, 'multiqc_fastqc.txt'),
        log:
            os.path.join(OD_LOG, 'multiqc_ref.log')
        threads: 1
        resources:
            mem_mb=20000
        params:
            title = 'Report using RefSeq annotations',
            comment = 'MultiQC report includes only RefSeq gene annotations.',
        singularity:
            config['MULTIQC_IMAGE']
        shell:
            """
            #unset LD_PRELOAD
            LC_ALL='en_US.utf-8'
            multiqc --no-megaqc-upload --verbose --force \
                --title '{params.title}' \
                --comment '{params.comment}' \
                --outdir {OD_MULTIQC} \
                --filename multiqc_report.html \
                {OD}/cutadapt \
                {OD}/fastqc \
                {OD}/log/*_Log.final.out \
                {OD}/metrics/*.refseq.RNAmetrics.txt \
                {OD}/fc/*.refseq.*  2> {log} \
            && mv -f {OD_MULTIQC}/multiqc_report_data/* {OD_MULTIQC}/ \
            && mv -f {OD_MULTIQC}/multiqc_report.html {output[0]} \
            && rmdir {OD_MULTIQC}/multiqc_report_data
            """
            
    rule multiqc_ens: # Rule for Ensembl annotations
        input:
            expand(os.path.join(OD_METRICS, '{sample}.ensembl.RNAmetrics.txt'), sample=sample_ids),
            expand(os.path.join(OD_FASTQC, '{name}_fastqc.html'), name=fastq_names_noext),
            expand(os.path.join(OD_FASTQC, '{name}_fastqc'), name=fastq_names_noext),
            expand(os.path.join(OD_STATS, '{sample}.bam.flagstat'), sample=sample_ids),
            expand(os.path.join(OD_STATS, '{sample}.bam.stats'), sample=sample_ids),
            expand(os.path.join(OD_STATS, '{sample}.bam.stats2'), sample=sample_ids),
            expand(os.path.join(OD_LOG, '{sample}_Log.final.out'), sample=sample_ids),
            expand(os.path.join(OD_FC, '{sample}.ensembl.cnt.summary'), sample=sample_ids),
            expand(os.path.join(OD_FC, '{sample}.ensembl.cnt.gz'), sample=sample_ids),
        output:
            os.path.join(OD, 'multiqc_report_ensembl.html'),
            os.path.join(OD_MULTIQC_ENS, 'multiqc_picard_RnaSeqMetrics.txt'),
            os.path.join(OD_MULTIQC_ENS, 'multiqc_featureCounts.txt'),
            os.path.join(OD_MULTIQC_ENS, 'multiqc_star.txt'),
            os.path.join(OD_MULTIQC_ENS, 'multiqc_fastqc.txt'),
        log:
            os.path.join(OD_LOG, 'multiqc_ens.log')
        threads: 1
        resources:
            mem_mb=20000
        params:
            title = 'Report using Ensembl annotations',
            comment = 'MultiQC report includes only Ensembl gene annotations.'
        singularity:
            config['MULTIQC_IMAGE']
        shell:
            """
            #unset LD_PRELOAD
            LC_ALL='en_US.utf-8'
            multiqc --no-megaqc-upload --verbose --force \
                --title '{params.title}' \
                --comment '{params.comment}' \
                --outdir {OD_MULTIQC_ENS} \
                --filename multiqc_report_ensembl.html \
                {OD}/fastqc \
                {OD}/cutadapt \
                {OD}/log/*_Log.final.out \
                {OD}/metrics/*.ensembl.RNAmetrics.txt \
                {OD}/fc/*.ensembl.*  2> {log} \
            && mv -f {OD_MULTIQC_ENS}/multiqc_report_ensembl_data/* {OD_MULTIQC_ENS}/ \
            && mv -f {OD_MULTIQC_ENS}/multiqc_report_ensembl.html {output[0]} \
            && rmdir {OD_MULTIQC_ENS}/multiqc_report_ensembl_data
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
