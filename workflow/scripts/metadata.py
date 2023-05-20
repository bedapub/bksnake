"""
Source:
    https://code.roche.com/PMDA/create_samples_table/-/blob/main/create_samples_table.py
"""
import sys
import pandas as pd
import os.path
from uuid import UUID
from pathlib import Path
import requests
import configparser

from .funcs import *


# Column names in metadata panda dataframe
ORGANISM_PROPERTY = 'Organism'
STRANDED_PROPERTY = 'Is Stranded (DAT)'
FIRST_READ_PROPERTY = 'FASTQ1'   # originally: FASTQ Forward (RDT)
SECOND_READ_PROPERTY = 'FASTQ2'  # originally: FASTQ Reverse (RDT)




"""
Login to MOOSE sample metadata database
"""
def login_and_get_headers_api(login_url, email, password, verify=True):
    login_data = {'email': email, 'password': password}
    login_res = requests.post(url=login_url, json=login_data, verify=verify)

    if login_res.status_code != 200:
        raise Exception(f'Login to API failed. {login_res.json()}')

    header = login_res.json()
    return header



"""
Get sample metadata from MOOSE database via API call
"""
def request_samples(api_login, api_password, dataset_id, study_id, api_url,
                    etl_api_url, verify=True, debug=False):
    
    if study_id is None or isBlank(study_id):
        endpoint = f'{etl_api_url}/api/studies/datasets/{dataset_id}/rdt/export'
    else:
        endpoint = f'{etl_api_url}/api/studies/id/{study_id}/datasets/id/{dataset_id}/rdt/export'
 
    params = {'format': 'tsv', 'use_cv_labels': 'true', 'prettify_headers': 'true'}
    if debug:
        print(f'Using endpoint {endpoint} and params {params}', file=sys.stdout)
    headers = login_and_get_headers_api( f'{api_url}/users/login', api_login, api_password, verify=verify )
    response = requests.get( endpoint, params=params, stream=True, headers=headers, verify=verify )
    if not response.ok:
        print(f'Server response:\n{response.text}', file=sys.stdout)
    response.raise_for_status()
    return response


"""
Read a configuration file for the MOOSE database
"""
def read_config(filename=None, debug=False):
    if filename is None:
        filename = Path.home() / CONFIG_PATH
    else:
        filename = Path(filename)
    if filename.is_file():
        if debug:
            print(f'Config file found: {filename}', file=sys.stderr)
        try:
            config = configparser.ConfigParser()
            config.read(filename)
            if 'default' not in config:
                print(f'Cannot find default section in the config file', file=sys.stderr)
            cfg_dict = dict(config['default'])
            if debug:
                print(f'Got config: {cfg_dict}', file=sys.stderr)
            return cfg_dict
        except Exception as e:
            print(f'Error loading config file {filename}: {e}', file=sys.stderr)

    return {}


"""
Updates:
    - if metadata_file exists, then read from that file.
    - allow combining multiple columns for group
    - remove special characters from group, e.g. greek symbols
    - remove first character from fastq files names if it is a slash '/'
    
    api_url: 'https://moose.roche.com/api'                                             
    etl_api_url: 'https://moose.roche.com/etl-api'      
"""
def get_metadata(dataset_id, study_id, group, credentials_file, api_url, etl_api_url, metadata_file):

    if os.path.exists(metadata_file):
        df = pd.read_csv(metadata_file, sep='\t')
    else:
        config = read_config(filename = credentials_file)
        request = request_samples(
                config.get('api_login'),
                config.get('api_password'),
                dataset_id,
                study_id,
                api_url,
                etl_api_url,
                verify = True)
        #df = pd.read_csv(request.raw, sep='\t')
        df = pd.read_csv(request.raw, sep='\t', dtype={'Readout ID (RDT)': 'object'})
 

    # check for duplicated columns
    if df.columns.duplicated().any():
        print( 'Warning: Got duplicated columns:',
                df.columns[df.columns.duplicated()].unique(),
  
              file=sys.stderr )
        
    # remove all NaN columns but keep 'Is Stranded (DAT)' column = STRANDED_PROPERTY
    col_names = df.columns
    df_copy = df
    df.dropna(axis='columns', how='all', inplace=True)
    if STRANDED_PROPERTY in df_copy:
        df[STRANDED_PROPERTY] = df_copy[STRANDED_PROPERTY]
    del df_copy

    
    # determine whether single- or paired-end reads
    if 'FASTQ Reverse (RDT)' in df:
        single_end = False
    else:
        single_end = True

        
    # rename some columns
    if single_end == True:
            df.rename(
                    columns={
                        'Readout ID (RDT)': '#ID',
                        'FASTQ Forward (RDT)': 'FASTQ1',
                        },
                    inplace=True)
    else:
        df.rename(
                columns={
                    'Readout ID (RDT)': '#ID',
                    'FASTQ Forward (RDT)': 'FASTQ1',
                    'FASTQ Reverse (RDT)': 'FASTQ2',
                    },
                inplace=True)

        
    # remove '/' if it is first character in fastq files names
    if 'FASTQ1' in df.columns:
        df['FASTQ1'] = df['FASTQ1'].str.replace('^/', '', regex=True)
    if 'FASTQ2' in df.columns:
        df['FASTQ2'] = df['FASTQ2'].str.replace('^/', '', regex=True)

        
    # find columns which are unique without parenthesis
    dup_cols = df.columns.str.replace(r' \([^(]+\)$', '', regex=True)
    dup_cols = dup_cols.duplicated(keep=False)
    
    
    # remove parenthesis suffixes among unique columns
    cols_to_rename = df.columns[~dup_cols]
    rename_dict = dict(
            zip(cols_to_rename, cols_to_rename.str.replace(r' \([^(]+\)$', '', regex=True))
            )
    df.rename(columns=rename_dict, inplace=True)

    
    # Determine the "group" name, may combine several columns
    if group == 'NA' or group == 'na' or group == 'n/a':
        assert 'GROUP' not in df.columns
        df['GROUP'] = 'noGroup' #'NA'
    elif isNotBlank(group):
        groupList = [clean_string(item) for item in group.split(';')]
        for col in groupList:
            df[col] = df[col].astype(str)
        try:
            df['GROUP'] = df[groupList].T.agg('_'.join)
            
            # Remove unicode characters in the GROUP column
            df['GROUP'] = df['GROUP'].str.encode('ascii', 'ignore').str.decode('ascii')
 
        except KeyError:
            print(f"Unknown column(s) '{group}'", file=sys.stderr)
            sys.exit(2)
    else:
        assert 'GROUP' not in df.columns
        df['GROUP'] = 'noGroup' #'NA'

    # Replace NA group name by 'noGroup'
    df['GROUP'] = df['GROUP'].replace(to_replace=r'^NaN$', value='noGroup', regex=True)
    df['GROUP'] = df['GROUP'].replace(to_replace=r"^nan$", value='noGroup', regex=True)
    df['GROUP'] = df['GROUP'].replace(to_replace=r"^NA$", value='noGroup', regex=True)
 
        
    # it may be that the ID are numbers, panda converts them to integers
    # but this is a problem for the subsequent statements,
    # thus we convert the whole column to string
    df['#ID'] = df['#ID'].astype(str)

    
    # escape spaces in #ID
    df['#ID'] = df['#ID'].str.replace(' ', '_', regex=True)

    if single_end == True:
        df = move_columns_left(df, ['#ID', 'GROUP', 'FASTQ1'])
    else:
        df = move_columns_left(df, ['#ID', 'GROUP', 'FASTQ1', 'FASTQ2'])


    # replace tabs and newlines by space
    df.replace('\t|\n|\r', ' ', regex=True, inplace=True)


    # replace white spaces and slashes in GROUP by underscore
    df['GROUP'] = df['GROUP'].str.replace(' ', '_', regex=False).replace('/', '_', regex=False)
    #df['GROUP'] = df['GROUP'].str.replace(' ', '_', regex=True).replace('/', '_', regex=True)

    return df




"""
Determine sequnecing library type: single-end or paired-end
"""
def library_type(meta):
    if FIRST_READ_PROPERTY in meta.columns and SECOND_READ_PROPERTY in meta.columns:
        return 'paired-end'
    elif FIRST_READ_PROPERTY in meta.columns:
        return 'single-end'
    else:
        raise Exception('\"'+FIRST_READ_PROPERTY+'\" property not present in metadata, \
                         thus cannot determine whether paired- or single-end sequencing library.')

        
"""
Determine strand orientation of the sequencing library

Normally this value should come from the Moose metadata as True or False

However, the strand may have 3 options:
 1. SECOND_READ_TRANSCRIPTION_STRAND (default, most common) 
 2. FIRST_READ_TRANSCRIPTION_STRAND 
 3. NONE
 
Strategy:
 1. if present and valid in config then use it from config
 2. if present and valid from moose/metadata then it from moose
 3. raise error

What if Moose/metadata 'Is Stranded (DAT)' is True but it is not 2nd but 1st read strand ? --> use config.yaml
"""
def strandedness(meta, config):

    default_strand = tuple(('SECOND_READ_TRANSCRIPTION_STRAND', 2, 2)) # most common option for paired-end library
    first_read = tuple(('FIRST_READ_TRANSCRIPTION_STRAND', 1, 1)) # less frequnent for paired-end library
    no_strand = tuple(('NONE', 0, 0)) # no strand or single-end library

    if library_type(meta) == 'single-end':
        return no_strand

    if 'library' in config and 'strand' in config['library']:
        sys.stderr.write('Sequencing library properties used from input config file.\n')
        if config['library']['strand'] == 'SECOND_READ_TRANSCRIPTION_STRAND':
            return default_strand
        elif config['library']['strand'] == 'FIRST_READ_TRANSCRIPTION_STRAND':
            return first_strand
        elif config['library']['strand'] == 'NONE':
            return no_strand
        else:
            raise Exception('Invalid sequencing library property in input config file: '+
                            config['library']['strand'])

    if STRANDED_PROPERTY in meta.columns:
        if len(meta[STRANDED_PROPERTY].unique()) == 1:
            if meta[STRANDED_PROPERTY].unique()[0]: # is stranded library
                return default_strand
            else:
                return no_strand
        else:
            raise Exception('Multiple values for sequencing library property \"'+
                            STRANDED_PROPERTY+'\" in metadata.\n')
    else:
        sys.stderr.write('Sequencing library property \"'+STRANDED_PROPERTY+
                         '\" not present in metadata nor in input config file. Thus, use \"unstranded\".\n')
        return no_strand
        
    return default_strand


"""
Determine fastq file locators from metadata
This depends whether paired- or single-end reads
"""
def fastq_locators(meta, config):
    
    if 'Raw' in meta:
        raw = meta['Raw'].str.strip('/').str.split('/', expand=True)
        nraw = raw.shape[1] - 1
    else:
        raise Exception('\"Raw\" column is missing in metadata.')
    
    if library_type(meta) == 'paired-end':
        data = {'file': meta['FASTQ1'],
                'uuid': raw[nraw],
                'name': meta['#ID']+'_1.fastq.gz'}
        df1 = pd.DataFrame(data)
        data = {'file': meta['FASTQ2'],
                'uuid': raw[nraw],
                'name': meta['#ID']+'_2.fastq.gz'}
        df2 = pd.DataFrame(data)
        locators = df1.append(df2)
        del [[df1],[df2]]
    else:
        data = {'file': meta['FASTQ1'],
                'uuid': raw[nraw],
                'name': meta['#ID']+'.fastq.gz'}
        df1 = pd.DataFrame(data)
        locators = df1
        del[df1]
    locators['locator'] = locators['uuid']+'/'+locators['file']
    locators.set_index('file', inplace=True)

    return locators


"""
Determine/update library type: paired-end or single-end
in the past this was given in the input config file
here, we overwrite the value in the config if present
by the value obtained from the metadata (via Fastq files)
"""
def update_library_type(config, meta):
    if 'library' in config:
        if 'type' in config['library']:
            config['library']['type'] = library_type(meta)
        else:
            config['library'] = ( {'type' : library_type(meta)} )
    else:
        config['library'] = ( {'type' : library_type(meta)} )
    return config


"""
Determine organism and update config dictionary

5  values from config to be updated:

1. species: 'hg38'
2. species_name: 'human'
3. organism: 'Homo sapiens'
4. genome_dir: '/path/to/genome/directory/hg38'
5. biokitr genes

What if the user wants a specific assembly, mm10 instead of mm39 ?
 if assembly version already specified in 'config' then do not overwrite.
"""
def update_organism(config, meta):
    
    if ORGANISM_PROPERTY in meta:
        if len(meta[ORGANISM_PROPERTY].unique()) == 1:
            org = map_organism_to_config(meta[ORGANISM_PROPERTY].unique()[0])
            if not 'species' in config:
                config['species'] = config['genomes'][org]['species']
            if not 'species_name' in config:
                config['species_name'] = config['genomes'][org]['species_name']
            if not 'organism' in config:
                config['organism'] = config['genomes'][org]['organism']
            if not 'genome_dir' in config:
                config['genome_dir'] = config['genomes'][org]['genome_dir']
            if 'biokitr' in config:
                if not 'genes' in config['biokitr']:
                    config['biokitr']['genes'] = config['genomes'][org]['biokitr_genes']
        else:
            raise Exception('Multiple values for \"'+
                            ORGANISM_PROPERTY+'\" property in metadata.\n')
        
    else:
        raise Exception('\"'+ORGANISM_PROPERTY+'\" property not present metadata.\n')
    
    return config


"""
Map Moose properties to organism/genome names in the common.yaml config file.
"""
def map_organism_to_config(organism):
    org = organism.lower()
    if org == 'human' or org == 'homo sapiens':
        return 'hg38'
    elif org == 'mouse' or org == 'mus musculus':
        return 'mm39'
    elif org == 'rat' or org == 'rattus norvegicus':
        return 'rn7'
    elif org == 'cyno' or org == 'cynomolgus' or org == 'macaca fascicularis' :
        return 'mfa5'
    elif org == 'pig' or org == 'sus scrofa':
        return 'ss11'
    elif org == 'rabbit' or org == 'oryctolagus cuniculus':
        return 'oc2'
    else:
        raise Exception('Unknown organism: \"'+organism+'\"')
    return org
