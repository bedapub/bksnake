"""
Calculate gene counts with feature counts from subreads
"""

"""
NOTE: Introduction of new parameter in v2.0.2 (https://subread.sourceforge.net): 
  
  --countReadPairs

New parameter '--countReadPairs' is added to featureCounts to explicitly specify that read pairs will be counted, 
and the '-p' option in featureCounts now only specifies if the input reads are paired end 
(it also implied that counting of read pairs would be performed in previous versions).
"""

# ------------------------------------------------------------------------------
rule fc:
    input:
        str = rules.strandness.output.txt,
        bam = rules.star.output.bam,
        gtf = rules.annotations.output.gtf,
    output:
        cnt = os.path.join(OD_FC,'{sample}.{db}.cnt.gz'),
        summary = os.path.join(OD_FC,'{sample}.{db}.cnt.summary')
    log:
        os.path.join(OD_LOG,'{sample}.fc_{db}.log')
    params:
        cnt = os.path.join(OD_FC,'{sample}.{db}.cnt'),
        opt = config['featurecounts']['parameters'],
    threads: 4
    resources:
        mem_mb = 30000
    singularity:
        config['SUBREAD_IMAGE']
    shell:
        """
        set -e

        paired=$(grep -c 'PairEnd Data' {input.str}) || true
        str=$(grep 'featurecounts=' {input.str} | sort | uniq | cut -d= -f2) || true

        if [[ "$paired" == 1 ]]; then
            p_param="-p --countReadPairs"
        else
            p_param=""
        fi

        # Run featureCounts and capture its exit status
        featureCounts \
            {params.opt} \
            ${{p_param}} \
            -s ${{str}} -a {input.gtf} \
            -T {threads} -o {params.cnt} {input.bam} 2> {log} || fc_failed=true

        # Handle the output files
        if [ -e {params.cnt} ]; then
            gzip -f {params.cnt}
            # Make sure the summary exists, or create it if missing
            if [ ! -e {output.summary} ]; then
                echo "Warning: Summary file not created by featureCounts. Creating empty file." >> {log}
                touch {output.summary}
            fi
        else
            echo "Warning: featureCounts did not produce expected output file. Creating empty files." >> {log}
            touch {output.cnt}
            touch {output.summary}
        fi

        # If featureCounts failed but we created the required output files,
        # the rule will succeed but we'll leave a note in the log
        if [ "${{fc_failed:-false}}" = true ]; then
            echo "ERROR: featureCounts failed but pipeline continued with empty output files." >> {log}
        fi
        """
