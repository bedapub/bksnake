"""
rules for getting raw fastq files from a local path
"""

# ------------------------------------------------------------------------------
# copy raw data from a local folder
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
#       loc = lambda wildcards, output: locators.loc[locators.index[locators['name'] == output[0].split('/')[2]].to_list()[0]]['locator']     
    threads: 1
    resources:
        mem_mb=1000
    shell:
        """
        cp -pr {params.loc} {output}
        """