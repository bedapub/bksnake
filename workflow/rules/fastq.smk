"""
rules for getting raw fastq files from a local path
"""

# ------------------------------------------------------------------------------
# make symbolic link of raw data to local folder
# the variable 'locators' is declared in Snakefile
rule get_fastq:
    input:
        os.path.join(OD_TMP, '{name}.input')
    output:
        os.path.join(OD_FASTQ, '{name}')
    log:
        os.path.join(OD_LOG, '{name}.get_fastq.log')
    params:
       loc = lambda wildcards, output: locators.loc[locators.index[locators['name'] == os.path.basename(output[0])].to_list()[0]]['locator']
    threads: 1
    resources:
        mem_mb=1000
    shell:
        """
        absolute_path=$(realpath "{params.loc}")
        ln -s ${{absolute_path}} {output}
        """