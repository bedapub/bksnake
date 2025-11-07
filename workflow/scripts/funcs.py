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
        # Construct the full path for the first read file
        fastq1_combination = os.path.join(row[FASTQ_FOLDER_PROPERTY], row[FIRST_READ_PROPERTY])
        if not os.path.exists(fastq1_combination):
            raise ValueError(f"{FIRST_READ_PROPERTY} path does not exist: {fastq1_combination}")
        
        # Check paths based on LIBRARY_LAYOUT
        if row["LIBRARY_LAYOUT"] == "PAIRED":
            # Construct the full path for the second read file
            fastq2_combination = os.path.join(row[FASTQ_FOLDER_PROPERTY], row[SECOND_READ_PROPERTY])
            if not os.path.exists(fastq2_combination):
                raise ValueError(f"{SECOND_READ_PROPERTY} path does not exist: {fastq2_combination}")
        elif row["LIBRARY_LAYOUT"] == "SINGLE":
            # Ensure SECOND_READ_PROPERTY is empty or NA for SINGLE layout
            if row[SECOND_READ_PROPERTY] not in ["", "NA", None, np.nan]:
                raise ValueError(f"For SINGLE layout, {SECOND_READ_PROPERTY} must be empty or NA.")
                
"""--------------------------------------------------------------------
Check for duplicates in the FASTQ1/FASTQ2 columns, considering LIBRARY_LAYOUT.
"""
def check_fastq_columns(df):
    """Check for duplicates in the FASTQ1/FASTQ2 columns, considering LIBRARY_LAYOUT."""
    required_columns = [FIRST_READ_PROPERTY, FASTQ_FOLDER_PROPERTY, "LIBRARY_LAYOUT"]
    if not all(column in df.columns for column in required_columns):
        raise ValueError(f"Missing required columns: {required_columns}")

    combinations = set()

    for index, row in df.iterrows():
        # Construct the full path for the first FASTQ file
        fastq1_combination = os.path.join(row[FASTQ_FOLDER_PROPERTY], row[FIRST_READ_PROPERTY])
        
        if row["LIBRARY_LAYOUT"] == "SINGLE":
            # For SINGLE layouts, check just the FASTQ1 column
            if fastq1_combination in combinations:
                raise ValueError(f"Duplicate combination found in SINGLE layout: {fastq1_combination}")
            combinations.add(fastq1_combination)
        
        elif row["LIBRARY_LAYOUT"] == "PAIRED":
            # For PAIRED layouts, check both FASTQ1 and FASTQ2 columns
            fastq2_combination = os.path.join(row[FASTQ_FOLDER_PROPERTY], row[SECOND_READ_PROPERTY])
            if fastq1_combination == fastq2_combination:
                raise ValueError(f"Identical FASTQ1 and FASTQ2 combinations found in PAIRED layout: {fastq1_combination}")
            
            if fastq1_combination in combinations:
                raise ValueError(f"Duplicate FASTQ1 combination found in PAIRED layout: {fastq1_combination}")
            if fastq2_combination in combinations:
                raise ValueError(f"Duplicate FASTQ2 combination found in PAIRED layout: {fastq2_combination}")
            
            combinations.add(fastq1_combination)
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
def library_layout(df):
    """ Check if the SECOND_READ_PROPERTY column exists."""
    if SECOND_READ_PROPERTY not in df.columns:
        # If the column doesn't exist, add it and set all row values to an empty string
        df[SECOND_READ_PROPERTY] = ""
    
    # Create a new column "LIBRARY_LAYOUT"
    def determine_library_layout(row):
        if row[SECOND_READ_PROPERTY] in ["", None, "NA", np.nan]: # Check if the column is empty, None, or NA
            return "SINGLE"
        else:
            return "PAIRED"
    
    df["LIBRARY_LAYOUT"] = df.apply(determine_library_layout, axis=1)
    return df
    

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

    # Determine library layout (single or paired) for each sample NEW
    library_layout(df)

    # Check Raw column
    check_raw_directories(df)
 
    # Check if the paths in FASTQ1 and SECOND_READ_PROPERTY columns exist
    check_fastq_paths(df)

    # Perform the FASTQ and Source Folder checks
    check_fastq_columns(df)
    # For debugging
    #print("DATAFRAME WITH LIBRARY LAYOUT IS:", file=sys.stderr)
    #print(df, file=sys.stderr)
    #raise ValueError('*** STOP PROGRAM FOR DEBUGGING ***')
    
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
Determine fastq file locators from metadata
"""
def fastq_locators(meta, config):
    """Determine fastq file locators from metadata using the LIBRARY_LAYOUT column.

    Parameters
    ----------
        meta : DataFrame
            Input DataFrame containing sample metadata, including the LIBRARY_LAYOUT column.

    Returns
    -------
        locators : DataFrame
            Fastq file locators derived from input metadata.    
    """
    # Check if required columns are present
    required_columns = [FASTQ_FOLDER_PROPERTY, FIRST_READ_PROPERTY, 'LIBRARY_LAYOUT']
    for col in required_columns:
        if col not in meta:
            raise ValueError(f'{col} column is missing in metadata.')

    # Ensure LIBRARY_LAYOUT column contains valid values
    if not meta['LIBRARY_LAYOUT'].isin(['SINGLE', 'PAIRED']).all():
        raise ValueError("LIBRARY_LAYOUT column must contain only 'SINGLE' or 'PAIRED' values.")

    # Split the FASTQ_FOLDER_PROPERTY column to extract the last component (uuid)
    raw = meta[FASTQ_FOLDER_PROPERTY].str.strip('/').str.split('/', expand=True)
    nraw = raw.shape[1] - 1  # Index of the last component of the split path

    # Initialize an empty list to store rows for the locators DataFrame
    rows = []

    # Iterate over each row in the metadata DataFrame
    for _, row in meta.iterrows():
        if row['LIBRARY_LAYOUT'] == 'PAIRED':
            # If the library is paired-end, create a row for each read (FASTQ1 and FASTQ2)
            rows.append({
                'file': row[FIRST_READ_PROPERTY],
                'uuid': raw.loc[_, nraw],
                'path': row[FASTQ_FOLDER_PROPERTY],
                'name': row['#ID'] + '_1.fastq.gz',
                'locator_uuid': f"{raw.loc[_, nraw]}/{row[FIRST_READ_PROPERTY]}",
                'locator': f"{row[FASTQ_FOLDER_PROPERTY]}/{row[FIRST_READ_PROPERTY]}"
            })
            rows.append({
                'file': row[SECOND_READ_PROPERTY],
                'uuid': raw.loc[_, nraw],
                'path': row[FASTQ_FOLDER_PROPERTY],
                'name': row['#ID'] + '_2.fastq.gz',
                'locator_uuid': f"{raw.loc[_, nraw]}/{row[SECOND_READ_PROPERTY]}",
                'locator': f"{row[FASTQ_FOLDER_PROPERTY]}/{row[SECOND_READ_PROPERTY]}"
            })
        elif row['LIBRARY_LAYOUT'] == 'SINGLE':
            # If the library is single-end, create a row only for the first read (FASTQ1)
            rows.append({
                'file': row[FIRST_READ_PROPERTY],
                'uuid': raw.loc[_, nraw],
                'path': row[FASTQ_FOLDER_PROPERTY],
                'name': row['#ID'] + '.fastq.gz',
                'locator_uuid': f"{raw.loc[_, nraw]}/{row[FIRST_READ_PROPERTY]}",
                'locator': f"{row[FASTQ_FOLDER_PROPERTY]}/{row[FIRST_READ_PROPERTY]}"
            })

    # Create the locators DataFrame from the list of rows
    locators = pd.DataFrame(rows)

    # Set the 'name' column as the index
    locators['name_as_index'] = locators['name']
    locators.set_index('name_as_index', inplace=True)

    return locators


"""------------------------------------------------------------------------------
"""
def update_config_with_genome_id(config, genome_id):
    """
    Update the config dictionary with the given genome ID.

    Parameters
    ----------
    config : dict
        Snakemake workflow configuration dictionary.
    genome_id : str
        The genome ID to update in the config.
    """
    genome_info = config['genomes'][genome_id]
    genome_info['genome_dir'] = os.path.normpath(os.path.join(config['genome_dir'], genome_info['genome_subdir']))
    config.setdefault('species', genome_info['species'])
    config.setdefault('species_name', genome_info['species_name'])
    config.setdefault('organism', genome_info['organism'])
    if 'biokitr' in config:
        config['biokitr'].setdefault('genes', genome_info['biokitr_genes'])

    
"""------------------------------------------------------------------------------
Several cases:


'organism' is a parameter coming from the metadata file
'species' is a parameter coming from the configuration file
'overrule_organism' is a new parameter coming from the configuration file
   the default is False. By this parameter, if organism does not match species, 
   then the species will be taken . For example, when mapping monkey data to 
   the human genome.

case 1: species and organism absent --> STOP
case 2: species only                --> OK
case 3: organism only               --> OK but choose default species
case 4: species and organism
        subcase 1: both match       --> OK
        subcase 2: no match and overrule is False --> STOP
        subcase 3: no mach and overrule is True   --> OK choose species

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
    # Case 1:
    if (GENOMEID_PROPERTY not in config or not config[GENOMEID_PROPERTY]) and ORGANISM_PROPERTY not in meta:
        raise ValueError(f"'{GENOMEID_PROPERTY}' property from config file and '{ORGANISM_PROPERTY}' property from metadata are invalid or missing.")

    # Case 2:
    elif GENOMEID_PROPERTY in config and ORGANISM_PROPERTY not in meta:
        genomeid = config.get(GENOMEID_PROPERTY)
        if genomeid not in config['genomes']:
            raise ValueError(f'The {GENOMEID_PROPERTY} parameter given in the config file {genomeid} is not in the genomes dictionary of the config file.')
        print(f"'{GENOMEID_PROPERTY}' property from config file but '{ORGANISM_PROPERTY}' property not in metadata.", file=sys.stderr)

    # Case 3:
    elif GENOMEID_PROPERTY not in config and ORGANISM_PROPERTY in meta:
        organism_from_metadata = check_uniqueness_of_organism(config, meta)  # otherwise raises an error!
        if organism_from_metadata not in [details['organism'] for details in config['genomes'].values()]:
            raise ValueError(f"'{ORGANISM_PROPERTY}' property from the metadata: {organism_from_metadata}, is not in config file.")
        genomeid = get_default_species_for_organism(config['genomes'], organism_from_metadata)
        print(f"'{ORGANISM_PROPERTY}' property from metadata: {organism_from_metadata}, but '{GENOMEID_PROPERTY}' property not given in config file. Use default species: {genomeid}.", file=sys.stderr)

    # Case 4:
    else:
        genomeid_from_config = config.get(GENOMEID_PROPERTY)
        organism_from_metadata = check_uniqueness_of_organism(config, meta)  # otherwise raises an error!
        organism_from_config = config['genomes'][genomeid_from_config]['organism']

        # Subcase 1:
        if organism_from_metadata == organism_from_config:
            genomeid = genomeid_from_config
            print(f"Organism from config and metadata match each other.", file=sys.stderr)

        # Subcase 2:
        elif 'overrule_organism' not in config or not config['overrule_organism']:
            raise ValueError(f"Organism from config: {organism_from_config} is different than organism from metadata file: {organism_from_metadata}. Use parameter 'overrule_organism' to override.")

        # Subcase 3:
        else:
            genomeid = genomeid_from_config
            print(f"Organism from config: {organism_from_config} overrules organism from metadata: {organism_from_metadata}.", file=sys.stderr)

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


"""------------------------------------------------------------------------------
"""
def check_single_default_per_organism(genomes):
    """
    Check if there is only one default: True value per organism in the genomes dictionary.

    This function iterates through the provided genomes dictionary and ensures that each organism
    has at most one genome marked as default (i.e., 'default': True). If more than one default genome
    is found for any organism, a ValueError is raised with details of the conflicting entries.

    Parameters:
    genomes (dict): Dictionary containing genome information. The dictionary should have the following structure:
                    {
                        'genome_id': {
                            'species': str,
                            'species_name': str,
                            'organism': str,
                            'genome_subdir': str,
                            'biokitr_genes': list,
                            'db': list,
                            'default': bool
                        },
                        ...
                    }

    Raises:
    ValueError: If more than one default: True value is found for any organism, with details of the conflicting entries.
    """
    organism_defaults = {}

    for genome, details in genomes.items():
        organism = details['organism']
        is_default = details.get('default', False)

        if organism not in organism_defaults:
            organism_defaults[organism] = []

        organism_defaults[organism].append((genome, is_default))

    for organism, entries in organism_defaults.items():
        default_count = sum(1 for _, is_default in entries if is_default)
        if default_count > 1:
            error_message = f"More than one default: True value found for organism '{organism}':\n"
            for genome, is_default in entries:
                error_message += f"  - Genome: {genome}, Default: {is_default}\n"
            raise ValueError(error_message)


"""------------------------------------------------------------------------------
"""            
def get_default_species_for_organism(genomes, organism):
    """
    Get the default species for a given organism from the genomes dictionary.

    This function iterates through the provided genomes dictionary and returns the species
    that is marked as default (i.e., 'default': True) for the specified organism. If no default
    species is found for the organism, or if the organism is not present in the dictionary,
    the function returns None.

    Parameters:
    genomes (dict): Dictionary containing genome information. The dictionary should have the following structure:
                    {
                        'genome_id': {
                            'species': str,
                            'species_name': str,
                            'organism': str,
                            'genome_subdir': str,
                            'biokitr_genes': list,
                            'db': list,
                            'default': bool
                        },
                        ...
                    }
    organism (str): The organism for which to find the default species.

    Returns:
    str or None: The default species for the given organism, or None if no default species is found.
    """
    for genome, details in genomes.items():
        if details['organism'] == organism and details.get('default', False):
            return details['species']
    return None



"""------------------------------------------------------------------------------
"""            
def check_uniqueness_of_organism(config, meta):
    """
    Check whether Organism values are unique in metadata
   
    meta : DataFrame
            Input DataFrame containing sample metadata.
    config : Dict
            Snakemake workflow configuration dictionary.

    Raises:
    ValueError: If organism are not unique.
    
    Returns:
    str : The unique, translated organism
    """
    unique_organisms = meta[ORGANISM_PROPERTY].unique()
    if len(unique_organisms) > 1:
        raise ValueError(f'The {ORGANISM_PROPERTY} column given in the metadata file contains multiple values.')
        
    organism_from_metadata = translate_species(unique_organisms[0], config['translations'], ignore_case=True)
        
    return organism_from_metadata
