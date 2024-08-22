import os
import sys
import csv
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
GENOMEID_PROPERTY = 'species'
ORGANISM_PROPERTY = 'Organism'
FIRST_READ_PROPERTY = 'FASTQ1'
SECOND_READ_PROPERTY = 'FASTQ2'
FASTQ_FOLDER_PROPERTY = 'Raw'

"""--------------------------------------------------------------------
default header used in the case when there is no header in the sample file.
"""
default_header = ('#ID', 'GROUP', FIRST_READ_PROPERTY, SECOND_READ_PROPERTY)


"""--------------------------------------------------------------------
Helper functions
"""
def clean_string(s):
    return s.strip()

def is_blank(my_string, check_not=False):
    return not (my_string and my_string.strip()) if not check_not else bool(my_string and my_string.strip())

def is_not_blank(myString):
    return bool(myString and myString.strip())

def move_columns_left(df, columns):
    cols = [col for col in columns if col in df.columns]
    return df[cols + [col for col in df.columns if col not in cols]]


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
    header = default_header
    samples = []
    with open(filename, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for elements in reader:
            if elements[0].startswith('#'):
                header = elements
            else:
                samples.append(OrderedDict(zip(header, elements)))
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
        writer = DictWriter(f, samples[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(samples)            
            

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
        sample_id = sample['#ID']
        if sample[SECOND_READ_PROPERTY] == 'none':
            fq = expand('{sample}', sample=sample_id)
            fastq_source[fq[0]] = os.path.abspath(sample[FIRST_READ_PROPERTY])
        else:
            fq = expand('{sample}_{end}', sample=sample_id, end=(1, 2))
            fastq_source[fq[0]] = os.path.abspath(sample[FIRST_READ_PROPERTY])
            fastq_source[fq[1]] = os.path.abspath(sample[SECOND_READ_PROPERTY])
        sample_fastq[sample_id] = fq
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
        #print('All sample ids are unique.')
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
    res = pd.read_table(filename, dtype={'#ID': str, 'GROUP': str})
    check_unique_sampleids(res.iloc[:, 0])

    res[FIRST_READ_PROPERTY] = res[FIRST_READ_PROPERTY].apply(
        lambda x: os.path.abspath(os.path.join(fileroot, x)) if not os.path.isabs(x) else x
    )

    res['GROUP'] = res['GROUP'].replace(to_replace=r'^(NaN|nan|NA)$', value='noGroup', regex=True)
    return res


"""------------------------------------------------------------------------------
"""
def biokit_sample_anno_to_phenoData(anno):
    """Convert biokit sample annotation DataFrame into the phenoData DataFrame
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
    anno['GROUP'] = anno['GROUP'].replace(to_replace=r'^(NaN|nan|NA)$', value='noGroup', regex=True)
    anno['ID_GROUP'] = anno['#ID'].astype(str) + "_" + anno['GROUP'].astype(str)
    anno = anno.drop(columns=[FIRST_READ_PROPERTY, SECOND_READ_PROPERTY, FASTQ_FOLDER_PROPERTY], errors='ignore')
    anno = anno.rename(columns={'#ID': 'ID'})
    columns_order = ['ID_GROUP', 'ID', 'GROUP'] + [col for col in anno.columns if col not in ['ID', 'GROUP', 'ID_GROUP']]
    return anno[columns_order]


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

    
"""--------------------------------------------------------------------
Check for paths of FASTQ1/2 files
"""
def check_fastq_paths(df):
    """Check if the paths in FIRST_READ_PROPERTY and SECOND_READ_PROPERTY columns exist."""
    for index, row in df.iterrows():
        fastq1_combination = os.path.join(row[FASTQ_FOLDER_PROPERTY], row[FIRST_READ_PROPERTY])
        if not os.path.exists(fastq1_combination):
            raise ValueError(f"{FIRST_READ_PROPERTY} path does not exist: {fastq1_combination}")

        if SECOND_READ_PROPERTY in df.columns:
            fastq2_combination = os.path.join(row[FASTQ_FOLDER_PROPERTY], row[SECOND_READ_PROPERTY])
            if not os.path.exists(fastq2_combination):
                raise ValueError(f"{SECOND_READ_PROPERTY} path does not exist: {fastq2_combination}")                

                
"""--------------------------------------------------------------------
Check for duplicates in the FASTQ1/2 columns
"""
def check_fastq_columns(df):
    required_columns = [FIRST_READ_PROPERTY, FASTQ_FOLDER_PROPERTY]
    if not all(column in df.columns for column in required_columns):
        raise ValueError(f"Missing required columns: {required_columns}")

    has_fastq2 = SECOND_READ_PROPERTY in df.columns
    combinations = set()

    for index, row in df.iterrows():
        fastq1_combination = os.path.join(row[FASTQ_FOLDER_PROPERTY], row[FIRST_READ_PROPERTY])
        if has_fastq2:
            fastq2_combination = os.path.join(row[FASTQ_FOLDER_PROPERTY], row[SECOND_READ_PROPERTY])
            if fastq1_combination == fastq2_combination:
                raise ValueError(f"Identical combinations found: {fastq1_combination}")

        if fastq1_combination in combinations:
            raise ValueError(f"Duplicate combination found: {fastq1_combination}")

        combinations.add(fastq1_combination)
        if has_fastq2:
            combinations.add(fastq2_combination)            
            

"""--------------------------------------------------------------------
"""
def check_raw_directories(df):           
    """Check if the directories in the FASTQ_FOLDER_PROPERTY column exist."""
    if FASTQ_FOLDER_PROPERTY not in df.columns:
        raise ValueError(f"Missing required column: {FASTQ_FOLDER_PROPERTY}")

    missing_dirs = [directory for directory in df[FASTQ_FOLDER_PROPERTY].unique() if not os.path.exists(directory)]
    if missing_dirs:
        raise ValueError(f"Directories do not exist: {', '.join(missing_dirs)}")            

    
"""--------------------------------------------------------------------
"""
def get_metadata_from_file(metadata_file, group):
    """Get metadata from user given flat file.

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
    if not os.path.exists(metadata_file):
        raise FileNotFoundError(f'Metadata file is missing. Abort! {metadata_file}')

    df = pd.read_csv(metadata_file, sep='\t')

    # Check Raw column
    check_raw_directories(df)

    # Check if the paths in FASTQ1 and SECOND_READ_PROPERTY columns exist
    check_fastq_paths(df)

    # Perform the FASTQ and Source Folder checks
    check_fastq_columns(df)

    # Check for duplicated columns
    if df.columns.duplicated().any():
        print('Warning: Got duplicated columns:', df.columns[df.columns.duplicated()].unique(), file=sys.stderr)

    # Remove all NaN columns
    df.dropna(axis='columns', how='all', inplace=True)

    # Determine whether single- or paired-end reads
    single_end = SECOND_READ_PROPERTY not in df

    # Remove '/' if it is first character in fastq files names
    df[FIRST_READ_PROPERTY] = df[FIRST_READ_PROPERTY].str.lstrip('/')
    if not single_end:
        df[SECOND_READ_PROPERTY] = df[SECOND_READ_PROPERTY].str.lstrip('/')

    # Find columns which are unique without parenthesis
    dup_cols = df.columns.str.replace(r' \([^(]+\)$', '', regex=True)
    dup_cols = dup_cols.duplicated(keep=False)

    # Remove parenthesis suffixes among unique columns
    cols_to_rename = df.columns[~dup_cols]
    rename_dict = dict(zip(cols_to_rename, cols_to_rename.str.replace(r' \([^(]+\)$', '', regex=True)))
    df.rename(columns=rename_dict, inplace=True)

    # Determine the "group" name, may combine several columns
    if is_blank(group) or group.lower() in ['na', 'n/a']:
        df['GROUP'] = 'noGroup'
    elif is_not_blank(group):
        if group not in df.columns:
            raise ValueError(f'The provided [metadata][group_name]="{group}" from the config file is not in the input metadata file: {metadata_file}. The columns in the metadata file are: {df.columns}')
        group_list = [clean_string(item) for item in group.split(';')]
        for col in group_list:
            df[col] = df[col].astype(str)
        try:
            df['GROUP'] = df[group_list].T.agg('_'.join)
            # Remove unicode characters in the GROUP column
            df['GROUP'] = df['GROUP'].str.encode('ascii', 'ignore').str.decode('ascii')
        except KeyError:
            print(f"Unknown column(s) '{group}'", file=sys.stderr)
            sys.exit(2)
    else:
        raise ValueError(f'GROUP error. Abort! [metadata][group_name]={group}')

    # Replace NA group name by 'noGroup'
    df['GROUP'] = df['GROUP'].replace(to_replace=r'^(NaN|nan|NA)$', value='noGroup', regex=True)

    # Convert ID column to string and escape spaces
    df['#ID'] = df['#ID'].astype(str).str.replace(' ', '_', regex=True)

    if single_end:
        df = move_columns_left(df, ['#ID', 'GROUP', FIRST_READ_PROPERTY])
    else:
        df = move_columns_left(df, ['#ID', 'GROUP', FIRST_READ_PROPERTY, SECOND_READ_PROPERTY])

    # Replace tabs and newlines by space
    df.replace('\t|\n|\r', ' ', regex=True, inplace=True)

    # Replace white spaces and slashes in GROUP by underscore
    df['GROUP'] = df['GROUP'].str.replace(' ', '_', regex=False).replace('/', '_', regex=False)

    # Add Gender column if missing
    normalized_columns = [col.lower() for col in df.columns]
    if 'sex' not in normalized_columns and 'gender' not in normalized_columns:
        print('NOTE: Columns "S/sex" or "G/gender" not present in sample metadata columns. Workflow will add column "Gender" with values "Unknown"', file=sys.stderr)
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
    if FASTQ_FOLDER_PROPERTY not in meta:
        raise ValueError(f'{FASTQ_FOLDER_PROPERTY} column is missing in metadata.')

    raw = meta[FASTQ_FOLDER_PROPERTY].str.strip('/').str.split('/', expand=True)
    nraw = raw.shape[1] - 1

    if library_type(meta) == 'paired-end':
        data1 = {
            'file': meta[FIRST_READ_PROPERTY],
            'uuid': raw[nraw],
            'path': meta[FASTQ_FOLDER_PROPERTY],
            'name': meta['#ID'] + '_1.fastq.gz'
        }
        data2 = {
            'file': meta[SECOND_READ_PROPERTY],
            'uuid': raw[nraw],
            'path': meta[FASTQ_FOLDER_PROPERTY],
            'name': meta['#ID'] + '_2.fastq.gz'
        }
        df1 = pd.DataFrame(data1)
        df2 = pd.DataFrame(data2)
        locators = pd.concat([df1, df2])
    else:
        data = {
            'file': meta[FIRST_READ_PROPERTY],
            'uuid': raw[nraw],
            'path': meta[FASTQ_FOLDER_PROPERTY],
            'name': meta['#ID'] + '.fastq.gz'
        }
        locators = pd.DataFrame(data)

    locators['locator_uuid'] = locators['uuid'] + '/' + locators['file']
    locators['locator'] = locators['path'] + '/' + locators['file']
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
    library_type_value = library_type(meta)
    config.setdefault('library', {}).setdefault('type', library_type_value)
    return config

"""------------------------------------------------------------------------------
"""
def update_config_with_genome_id(config, genome_id):
    """Update the config dictionary with genome information based on genome_id."""       
    genome_info = config['genomes'][genome_id]
    genome_info['genome_dir'] = os.path.normpath(os.path.join(config['genome_dir'], genome_info['genome_subdir']))
    config.setdefault('species', genome_info['species'])
    config.setdefault('species_name', genome_info['species_name'])
    config.setdefault('organism', genome_info['organism'])
    if 'biokitr' in config:
        config['biokitr'].setdefault('genes', genome_info['biokitr_genes'])

    
"""------------------------------------------------------------------------------
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
    genomeid_from_config = config.get(GENOMEID_PROPERTY)
    organism_from_config = None

    if genomeid_from_config:
        if genomeid_from_config not in config['genomes']:
            raise ValueError(f'The {GENOMEID_PROPERTY} parameter given in the config file {genomeid_from_config} is not in the genomes dictionary of the config file.')
        organism_from_config = config['genomes'][genomeid_from_config]['organism']

    organism_from_metadata = None
    genomeid_from_metadata = None

    if ORGANISM_PROPERTY in meta:
        unique_organisms = meta[ORGANISM_PROPERTY].unique()
        if len(unique_organisms) > 1:
            raise ValueError(f'The {ORGANISM_PROPERTY} column given in the metadata file contains multiple values.')

        organism_from_metadata = translate_species(unique_organisms[0], config['translations'], ignore_case=True)
        if not genomeid_from_config:
            check_organism_in_config(config, organism_from_metadata)

    if organism_from_config and organism_from_metadata and organism_from_config != organism_from_metadata:
        raise ValueError(f'Organism from config: {organism_from_config} is different than organism from metadata file: {organism_from_metadata}.')

    genomeid = genomeid_from_config or genomeid_from_metadata
    if not genomeid:
        raise ValueError(f'Genome id missing in config ({GENOMEID_PROPERTY} parameter) or metadata file ({ORGANISM_PROPERTY} column).')

    update_config_with_genome_id(config, genomeid)
    return config        


"""------------------------------------------------------------------------------
"""
def check_organism_in_config(config, organism):
    """Check whether the value of 'organism' is in the config dictionary and whether it is unique or occurs multiple times.

    Parameters
    ----------
    config : dict
        Configuration dictionary containing genomes data.
    organism : str
        The organism value to check in the config dictionary.

    Returns
    -------
    None
    """
    # Initialize a counter for the organism
    organism_count = 0

    # Iterate through the genomes to count occurrences of the organism
    for genomeid, genome_info in config['genomes'].items():
        if genome_info['organism'].lower() == organism.lower():
            organism_count += 1

    # Check if the organism is unique or occurs multiple times
    if organism_count == 0:
        raise Exception(f'The organism "{organism}" is not found in the config dictionary.')
    elif organism_count == 1:
        print(f'The organism "{organism}" is unique in the config dictionary.', file=sys.stderr)
    else:
        raise Exception(f'The organism "{organism}" occurs multiple times ({organism_count}) in the config dictionary. Use the config file parameter {GENOMEID_PROPERTY} to specify.')


"""------------------------------------------------------------------------------
"""
def get_translation_dict(config):
    """Get the translation dictionary from the config dictionary."""
    return config['translations']


"""------------------------------------------------------------------------------
"""
def translate_species(input_value, translation_dict, ignore_case=True):
    """Translate input species value to standardized output value, with optional case sensitivity.

    Parameters
    ----------
    input_value : str
        The input species value to translate.
    translation_dict : dict
        The translation dictionary loaded from the YAML file.
    ignore_case : bool, optional
        Whether to ignore case when matching (default is True).

    Returns
    -------
    str
        The standardized output value if found, otherwise the original input value.
    """
    if ignore_case:
        input_value_lower = input_value.lower()
        for standard_value, synonyms in translation_dict.items():
            if input_value_lower in (synonym.lower() for synonym in synonyms):
                return standard_value
    else:
        for standard_value, synonyms in translation_dict.items():
            if input_value in synonyms:
                return standard_value
    return input_value
