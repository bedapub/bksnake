import os
import sys
import gzip
import subprocess # to get git hash
import pandas as pd
import numpy as np

from os.path import dirname, join
from collections import OrderedDict
from csv import DictWriter
from snakemake.io import expand



"""--------------------------------------------------------------------
Column names in metadata panda dataframe
"""
ORGANISM_PROPERTY = 'Organism'
FIRST_READ_PROPERTY = 'FASTQ1'
SECOND_READ_PROPERTY = 'FASTQ2'

"""--------------------------------------------------------------------
default header used in the case when there is no header in the sample file.
"""
default_header = ('#ID', 'GROUP', 'FASTQ1', 'FASTQ2')


"""--------------------------------------------------------------------
Helper functions
"""
def clean_string(s):
    return s.strip()

def isBlank (s):
    if s and s.strip():
        return False
    return True

def isNotBlank (s):
    if s and s.strip():
        return True
    return False

def move_columns_left(df, columns):
    df_left = df[columns]
    df_right = df[[col for col in df.columns if col not in df_left.columns]]
    df = pd.concat([df_left, df_right], axis='columns')
    return df


"""--------------------------------------------------------------------
"""
def read_samples(filename):
    """Parses samples.txt into list of dictionaries for every sample.
    Parameters
    ----------
    filename : str
        Name of tab-separated file with a header.
    Returns
    -------
    list of collections.OrderedDict
        Every item in the represents a single line from the tab-separated file.
    """
    global default_header
    header = default_header

    samples = []
    for line in open(filename):
        elements = line.rstrip('\n').split('\t')
        if line.startswith('#'):
            header = elements
        else:
            samples.append(OrderedDict(zip(
                header,
                elements
            )))
    return samples


"""--------------------------------------------------------------------
"""
def write_samples_file(samples, filename):
    """Recreates samples.txt from a list of ordered dictionaries.
    Parameters
    ----------
    samples : list of collections.OrderedDict
        every item is a sample.
    filename : str
        output filename to export tab-separated file.
    """
    assert len(samples) > 0
    with open(filename, 'w') as f:
        w = DictWriter(f, samples[0].keys(), delimiter='\t')
        w.writeheader()
        for sample in samples:
            w.writerow(sample)


"""--------------------------------------------------------------------
"""
def get_fastq_deps(samples):
    """Create dependency dictionaries for every sample.
    Parameters
    ----------
    samples : list of collections.OrderedDict
        Every element corresponds to a sample.
    Returns
    -------
    fastq_source : dict of str: str
        Mapping of fastq filenames to their paths. E.g.:
        # single end
        {"ABC123": "/data/samples/sample/ABC123.fastq.gz"}
        # pair end
        {
            "ABC123_1": "/data/samples/sample/ABC123_1.fastq.gz",
            "ABC123_2": "/data/samples/sample/ABC123_2.fastq.gz",
        }
    sample_fastq : dict of str: list of str
        Mapping of sample ids to fastq file names. E.g.:
        {"ABC123": ["ABC123_1", "ABC123_2"]}
    """
    fastq_source = {}
    sample_fastq = {}

    for sample in samples:
        if sample['FASTQ2'] == 'none':
            # single ended
            fq = expand('{sample}', sample=sample['#ID'])
            fastq_source[fq[0]] = os.path.abspath(sample['FASTQ1'])
        else:
            # paired end
            fq = expand('{sample}_{end}', sample=sample['#ID'], end=(1, 2))
            fastq_source[fq[0]] = os.path.abspath(sample['FASTQ1'])
            fastq_source[fq[1]] = os.path.abspath(sample['FASTQ2'])

        sample_fastq[sample['#ID']] = fq
    return fastq_source, sample_fastq



"""--------------------------------------------------------------------
Determine the git repository commit hash to create git url with hash
"""
def get_git_revision_short_hash() -> str:
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()




"""------------------------------------------------------------------------------
Check whether the input sample IDs are unique or not    
"""
def check_unique_sampleids(ids):
    """Check whether the input sample IDs are unique or not using pandas.

    Parameters
    ----------
    ids : list or array-like
        Input sample identifiers.
    """
    ids_series = pd.Series(ids)
    dups = ids_series[ids_series.duplicated(keep=False)]

    if not dups.empty:
        print('Duplicated ids: ' + ', '.join(dups.unique()), file=sys.stderr)
        print('Duplicated sample ids detected. Program quits.', file=sys.stderr)
        return False
    else:
        print('All sample ids are unique.')
        return True
    

"""------------------------------------------------------------------------------
Read sample annotations from a flat file
"""
def read_biokit_sample_anno(filename):
    """Read biokit sample annotation file, ensure that file
    paths are absolute paths, and returns a DataFrame
    
    Parameters
    ----------
    filename : str
        Input file name
        
    Returns
    -------
    res : DataFrame
        DataFrame containing sample information
    """
    fileroot = os.path.abspath(os.path.dirname(filename))
    res = pd.read_table(filename, dtype={'#ID':str, 'GROUP':str})
    check_unique_sampleids(res.iloc[:, 0])

    for filename in res.iloc[:,2]:
        if not os.path.isabs(filename):
            absfilename = os.path.abspath(os.path.join(fileroot, filename))
            res.iloc[:,2].replace({filename:absfilename}, inplace=True)
#    replace NA in GROUP column 
    res['GROUP'] = res['GROUP'].astype(str)
    res['GROUP'] = res['GROUP'].replace(to_replace=r'^NaN$', value='noGroup', regex=True)
    res['GROUP'] = res['GROUP'].replace(to_replace=r'^nan$', value='noGroup', regex=True)
    res['GROUP'] = res['GROUP'].replace(to_replace=r'^NA$', value='noGroup', regex=True)

    return(res)


"""------------------------------------------------------------------------------
Convert sample annotations to PhenoData format
"""
def biokit_sample_anno_to_phenoData(anno):
    """convert biokit sample annotation DataFrame into the phenoData DataFrame
    in the same format as the biokit output pipeline.

    Parameters
    ----------
        anno : DataFrame 
            Metadata from the function 'read_biokit_sample_anno'

    Returns
    -------
        res : DataFrame 
            A DataFrame object with the first three columns named ID_GROUP, ID, and
            GROUP. The rest columns contain sample annotation information in the
            input annotation file
    """
    anno['GROUP'] = anno['GROUP'].astype(str)
    anno['GROUP'] = anno['GROUP'].replace(to_replace=r'^NA$', value='noGroup', regex=True)

    id_group = anno.iloc[:,0].astype(str) + "_" + anno.iloc[:,1].astype(str)
    keepcols = [True] * len(anno.columns)
    keepcols[2:4] = [False]*2 # FASTQ1 and FASTQ2
    rest = anno.iloc[:, keepcols].copy()
    res = pd.concat([rest, id_group], axis=1)
    res.rename(columns={'#ID': 'ID', 0:'ID_GROUP'}, inplace=True)

    return(res)


"""------------------------------------------------------------------------------
Translate input biokit sample annotation file to output phenoData.meta
used, for example, in workflow/Snakefile
"""
def translate_biokit_to_phenoData_meta(infile, outfile):
    """Translate input biokit sample annotation file to output phenoData.meta
 
    Parameters
    ----------
        infile : str
            Input sample annotation file name.
        output : str
            Output phenoData.meta file name.
    """
    anno = read_biokit_sample_anno(infile)
    outdf = biokit_sample_anno_to_phenoData(anno)
    outdf.to_csv(outfile, sep='\t', index=False)


"""------------------------------------------------------------------------------
Get metadata from flat file or from study and sample metadata MongoDB database.
Updates:
    - if metadata_file exists, then read from that file.
    - allow combining multiple columns for group
    - remove special characters from group, e.g. greek symbols
    - remove first character from fastq files names if it is a slash '/'
    - add Gender column if missing
"""
def get_metadata_from_file(metadata_file, group):
    """Get metadata from flat file or from study and sample metadata MongoDB database.
    
    Parameters
    ----------
        metadata_file : str
            Input file name containing sample metadata information.
        group : str
            Name of the group or condition variable, ie one of the header/column names
            
    Returns
    -------
        df : DataFrame
            Modified and reformatted DataFrame containing sample metadata.
    """
    if os.path.exists(metadata_file):
        df = pd.read_csv(metadata_file, sep='\t')
    else:
        raise Exception('Metadata file is missing. Abort! '+metadata_file)

    # check for duplicated columns
    if df.columns.duplicated().any():
        print( 'Warning: Got duplicated columns:',
                df.columns[df.columns.duplicated()].unique(),  
              file=sys.stderr )

        
    # remove all NaN columns
    col_names = df.columns
    df_copy = df
    df.dropna(axis='columns', how='all', inplace=True)
    del df_copy


    # determine whether single- or paired-end reads
    if 'FASTQ2'in df:
        single_end = False
    else:
        single_end = True
 
        
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
    if isBlank(group) or group == 'NA' or group == 'na' or group == 'n/a':
        #assert 'GROUP' not in df.columns
        df['GROUP'] = 'noGroup'
    elif isNotBlank(group):
        if group not in df.columns:
            raise Exception('The provided [metadata][group_name]="'
                            +group
                            +'" from the config file is not in the input metadata file: '
                            +metadata_file+'. The colums are in the metadata file are: '
                            +str(df.columns))
        #assert 'GROUP' not in df.columns
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
        raise Exception('GROUP error. Abort! [metadata][group_name]='+group)
       

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
    
    # add Gender column if missing, could be used for biokitr package (in development)
    normalized_columns = [col.lower() for col in df.columns]
    if 'sex' not in normalized_columns and 'gender' not in normalized_columns:
        print(f'NOTE: Columns "S/sex" or "G/gender" not present in sample metadata columns. Add column "Gender" with values "Unknown"', file=sys.stderr)
        df['Gender'] = 'Unknown'

    return df



"""------------------------------------------------------------------------------
Determine sequnecing library type: single-end or paired-end
"""
def library_type(meta):
    """Determine sequnecing library type: single-end or paired-end.
    
    Parameters
    ----------
        meta : DataFrame
            Input DataFrame containing sample metadata.
    
    Returns
    -------
        library_type : str
            'paired-end' or 'single-end' derived from input metadata.
    """
    if FIRST_READ_PROPERTY in meta.columns and SECOND_READ_PROPERTY in meta.columns:
        return 'paired-end'
    elif FIRST_READ_PROPERTY in meta.columns:
        return 'single-end'
    else:
        raise Exception('\"'+FIRST_READ_PROPERTY+'\" property not present in metadata, \
                         thus cannot determine whether paired- or single-end sequencing library.')



"""------------------------------------------------------------------------------
Determine fastq file locators from metadata
This depends whether paired- or single-end reads
"""
def fastq_locators(meta, config):
    """Determine fastq file locators from metadata.
    
    Parameters
    ----------
        meta : DataFrame
            Input DataFrame containing sample metadata.
        config : Dict
            Snakemake workflow configuration dictionary.

    Returns
    -------
        locators : DataFrame
            Fastq file locators derived from input metadata.    
    """
    
    if 'Raw' in meta:
        raw = meta['Raw'].str.strip('/').str.split('/', expand=True)
        nraw = raw.shape[1] - 1
    else:
        raise Exception('\"Raw\" column is missing in metadata.')
    
    if library_type(meta) == 'paired-end':
        data = {'file': meta['FASTQ1'],
                'uuid': raw[nraw],
                'path': meta['Raw'],               
                'name': meta['#ID']+'_1.fastq.gz'}
        df1 = pd.DataFrame(data)
        data = {'file': meta['FASTQ2'],
                'uuid': raw[nraw],
                'path': meta['Raw'],
                'name': meta['#ID']+'_2.fastq.gz'}
        df2 = pd.DataFrame(data)
        locators = pd.concat([df1, df2])
        del df1, df2
    else:
        data = {'file': meta['FASTQ1'],
                'uuid': raw[nraw],
                'path': meta['Raw'],
                'name': meta['#ID']+'.fastq.gz'}
        df1 = pd.DataFrame(data)
        locators = df1
        del df1
    # old:
    #locators['locator'] = locators['uuid']+'/'+locators['file']
    #locators.set_index('file', inplace=True)
    # new:
    locators['locator_uuid'] = locators['uuid']+'/'+locators['file']
    locators['locator'] = locators['path']+'/'+locators['file']
    locators['name_as_index'] = locators['name']
    locators.set_index('name_as_index', inplace=True)

    return locators



"""------------------------------------------------------------------------------
"""
def update_library_type(config, meta):
    """Determine/update library type: paired-end or single-end in the past
    this was given in the input config file. Here, we overwrite the value 
    in the config if present by the value obtained from the metadata 
    (via Fastq files).
    
    Parameters
    ----------
        meta : DataFrame
            Input DataFrame containing sample metadata.
        config : Dict
            Snakemake workflow configuration dictionary.
    
    Returns
    -------
        config : Dict
            Updated Snakemake workflow configuration dictionary.
    """
    if 'library' in config:
        if 'type' in config['library']:
            config['library']['type'] = library_type(meta)
        else:
            config['library'] = ( {'type' : library_type(meta)} )
    else:
        config['library'] = ( {'type' : library_type(meta)} )
    return config



"""------------------------------------------------------------------------------
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
    """Determine organism and update config dictionary
    
    Parameters
    ----------
        meta : DataFrame
            Input DataFrame containing sample metadata.
        config : Dict
            Snakemake workflow configuration dictionary.
    
    Returns
    -------
        config : Dict
            Updated Snakemake workflow configuration dictionary.      
    """  
    if ORGANISM_PROPERTY in meta:
        if len(meta[ORGANISM_PROPERTY].unique()) == 1:
            org = map_organism_to_config(meta[ORGANISM_PROPERTY].unique()[0])
            config['genomes'][org]['genome_dir'] = os.path.normpath(os.path.join(config['genome_dir'], config['genomes'][org]['genome_subdir']))
            if not 'species' in config:
                config['species'] = config['genomes'][org]['species']
            if not 'species_name' in config:
                config['species_name'] = config['genomes'][org]['species_name']
            if not 'organism' in config:
                config['organism'] = config['genomes'][org]['organism']
            if 'biokitr' in config:
                if not 'genes' in config['biokitr']:
                    config['biokitr']['genes'] = config['genomes'][org]['biokitr_genes']
        else:
            raise Exception('Multiple values for \"'+
                            ORGANISM_PROPERTY+'\" property in metadata.\n')
    else:
        raise Exception('\"'+ORGANISM_PROPERTY+'\" property not present metadata.\n')
    
    return config



"""------------------------------------------------------------------------------
Map species properties to organism/genome names in the common.yaml config file.
"""
def map_organism_to_config(organism):
    """Map species properties to organism/genome names in the common.yaml config file.
    
    Parameters
    ----------
        organism : str
            Long name of the organism, e.g. 'human' or 'Homo sapiens'.
    
    Returns
    -------
        org : str
            Short name of the organism, e.g. 'hg38'.
    """
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
    elif org == 'green monkey' or org == 'african green monkey' or org == 'chlorocebus sabaeus':
        return 'Vero_WHO_p1.0'
    elif org == 'rhesus' or org == 'rhesus monkey' or org == 'macaca mulatta':
        return 'Mmul_10'
    elif org == 'hamster' or org == 'chinese hamster' or org == 'cricetulus griseus':
        return 'CriGri_PICRH_1_0'       
    elif org == 'dog' or org == 'canis lupus familiaris' or org == 'canis lupus':
        return 'ROS_Cfam_1.0 '       
    
# new: commented out in order to allow for 'custom' genomes
#    else:
#        raise Exception('Unknown organism: \"'+organism+'\"')
    return organism
