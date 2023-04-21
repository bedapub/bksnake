"""
rules for getting raw fastq files from arvados, ross or local path
"""

# Variables
if 'ross' in config:
    ross_profile = config['ross']['aws_profile']  # 'RSIROSS-BIM'
    ross_url = config['ross']['aws_endpoint_url'] # 'https://ross.science.roche.com'
    ross_bucket = config['ross']['aws_bucket']    # 'pred-bioinfo-bi-reference'
else:
    ross_profile = 'n/a'
    ross_url = 'n/a'
    ross_bucket = 'n/a'
location_raw = 'local'


# Process the metadata dataframe 'meta'

if 'ross' in config and ross_url+'/'+ross_bucket in list(set(meta['Raw']))[0]:
    location_raw = 'ross'
    if config['library']['type'] == 'paired-end':
        data = {'file': meta['FASTQ1'],
                's3uri': meta['Raw'].str.replace(ross_url,'s3:/'),
                'name': meta['#ID']+'_1.fastq.gz'}
        df1 = pd.DataFrame(data)
        data = {'file': meta['FASTQ2'],
                's3uri': meta['Raw'].str.replace(ross_url,'s3:/'),
                'name': meta['#ID']+'_2.fastq.gz'}
        df2 = pd.DataFrame(data)
        locators = df1.append(df2)
        del [[df1],[df2]]
    else:
        data = {'file': meta['FASTQ1'],
                's3uri': meta['Raw'].str.replace(ross_url,'s3:/'),
                'name': meta['#ID']+'.fastq.gz'}
        df1 = pd.DataFrame(data)
        locators = df1
        del[df1]
    locators['locator'] = locators['s3uri']+'/'+locators['file']
    locators.set_index('file', inplace=True)

elif ('https://wb2.arkau.roche.com' in list(set(meta['Raw']))[0]) | \
     ('https://wb2.arind.roche.com' in list(set(meta['Raw']))[0]) | \
     ('https://wb2.arsha.roche.com' in list(set(meta['Raw']))[0]):
    location_raw = 'arvados'
    if config['library']['type'] == 'paired-end':
        data = {'file': meta['FASTQ1'],
                'uuid': meta['Raw'].str.split('/', expand=True)[4],
                'name': meta['#ID']+'_1.fastq.gz'}
        df1 = pd.DataFrame(data)
        data = {'file': meta['FASTQ2'],
                'uuid': meta['Raw'].str.split('/', expand=True)[4],
                'name': meta['#ID']+'_2.fastq.gz'}
        df2 = pd.DataFrame(data)
        locators = df1.append(df2)
        del [[df1],[df2]]
    else:
        data = {'file': meta['FASTQ1'],
                'uuid': meta['Raw'].str.split('/', expand=True)[4],
                'name': meta['#ID']+'.fastq.gz'}
        df1 = pd.DataFrame(data)
        locators = df1
        del[df1]
    locators['locator'] = locators['uuid']+'/'+locators['file']
    locators.set_index('file', inplace=True)

else:
    location_raw = 'local'
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


# rules
# ---------------------------------------------------------------
if location_raw == 'arvados':
    rule get_fastq:
        input:
            os.path.join(OD_TMP, '{name}.input')
        output:
            os.path.join(OD_FASTQ, '{name}')
        log:
            os.path.join(OD_LOG, '{name}.get_fastq.log')
        params:
           loc = lambda wildcards, output: locators.loc[locators.index[locators['name'] == output[0].split('/')[2]].to_list()[0]]['locator'],
           uuid = lambda wildcards, output: locators.loc[locators.index[locators['name'] == output[0].split('/')[2]].to_list()[0]]['uuid'],
           file_list = os.path.join(OD_FASTQ, '{name}.arv_file_list.tmp')
        threads: 1
        resources:
            mem_mb=1000
        singularity:
            config['ARV_CLI_IMAGE']
        shell:
            """
            arv-ls {params.uuid} > {params.file_list}
            file=`echo {params.loc} | sed 's/{params.uuid}\///g'`
            n=`cat {params.file_list} | awk 'BEGIN{{n=0}}{{ if ($1 == "./'$file'") {{n=n+1}}}}END{{print n}}'`
            if [[ $n -eq 1 ]]; then
                arv-get --no-progress {params.loc} > {output} && chmod -R -x {output} && chmod -R gu+rw {output} && rm {params.file_list}
            else
                printf "\n\nError in rule get fastq files from Arvados: Fastq file {params.loc} is missing or has a wrong name! \n\
                Please check the readouts in the Moose data set and/or the fastq files names in the Arvados collection with UUID {params.uuid}.\n\
                Note that in Moose the 'FASTQ Forward' and 'FASTQ Reverse' fields must match the files names in the corresponding Arvados collection.\n\
                The required fastq file {params.loc} is not present in the Arvados collection with UUID {params.uuid}.\n\
                The following files are present in the Arvados collection:\n"
                cat {params.file_list}
            fi
            """
elif location_raw == 'ross':
    rule get_fastq:
        input:
            os.path.join(OD_TMP, '{name}.input')
        output:
            (os.path.join(OD_FASTQ, '{name}'))
        log:
            os.path.join(OD_LOG, '{name}.get_fastq.log')
        params:
           loc = lambda wildcards, output: locators.loc[locators.index[locators['name'] == output[0].split('/')[2]].to_list()[0]]['locator'],
           url = ross_url,
           profile = ross_profile
        threads: 1
        resources:
            mem_mb=1000
        singularity:
            config['AWS_CLI_IMAGE']
        shell:
            """
            aws s3 cp {params.loc} {output} --endpoint-url {params.url} --profile {params.profile}
            """
else: # copy from local folder
    rule get_fastq:
        input:
            os.path.join(OD_TMP, '{name}.input')
        output:
            (os.path.join(OD_FASTQ, '{name}'))
        log:
            os.path.join(OD_LOG, '{name}.get_fastq.log')
        params:
           loc = lambda wildcards, output: locators.loc[locators.index[locators['name'] == output[0].split('/')[2]].to_list()[0]]['locator'],
        threads: 1
        resources:
            mem_mb=1000
        shell:
            """
            cp -pr {params.loc} {output}
            """
