"""
rules for getting raw fastq files from a local path
"""

# ------------------------------------------------------------------------------
# Process the metadata dataframe 'meta'
if config['library']['type'] == 'paired-end':
    data = {'file': meta['FASTQ1'],
            'path': meta['Raw'],
            'name': meta['#ID']+'_1.fastq.gz'}
    df1 = pd.DataFrame(data)
    data = {'file': meta['FASTQ2'],
            'path': meta['Raw'],
            'name': meta['#ID']+'_2.fastq.gz'}
    df2 = pd.DataFrame(data)
    locators = df1.append(df2)
    del [[df1],[df2]]
else:
    data = {'file': meta['FASTQ1'],
            'path': meta['Raw'],
            'name': meta['#ID']+'.fastq.gz'}
    df1 = pd.DataFrame(data)
    locators = df1
    del[df1]
locators['locator'] = locators['path']+'/'+locators['file']
locators.set_index('file', inplace=True)

# ------------------------------------------------------------------------------
# copy raw data from a local folder
rule get_fastq:
    input:
        os.path.join(OD_TMP, '{name}.input')
    output:
        (os.path.join(OD_FASTQ, '{name}'))
    log:
        os.path.join(OD_LOG, '{name}.get_fastq.log')
    params:
       loc = lambda wildcards, output: locators.loc[locators.index[locators['name'] == output[0].split('/')[2]].to_list()[0]]['locator']
    threads: 1
    resources:
        mem_mb=1000
    shell:
        """
        cp -pr {params.loc} {output}
        """