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



# --------------------------------------------------------------------
default_header = ('#ID', 'GROUP', 'FASTQ1', 'FASTQ2')
"""str: default header used in the case when there is no header in the sample
file.
"""

# --------------------------------------------------------------------
# Helper functions
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


# --------------------------------------------------------------------
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


# --------------------------------------------------------------------
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


# --------------------------------------------------------------------
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


# --------------------------------------------------------------------
# Function to create a data table with mapping statistics from STAR
def star_mapping_stats (outfile, phenoFile, indir):
    print("star_mapping_stats")
    pheno = pd.read_csv(phenoFile, header = 0, sep = '\t', usecols=['ID','GROUP'], 
            dtype={'ID':str,'GROUP':str},
            keep_default_na=False, na_values=[''])

    print(pheno)
    # Read in all files with star output
    i = 0
    for sample in list(pheno['ID']):
        f = os.path.join(indir, sample+'_Log.final.out')
        df = pd.read_csv(f, header = None, sep = '\t')
#        s =  f.replace('_Log.final.out', '')
        df.columns = ['0', sample]
        if i == 0:
            tab = df
        else:
            tab = pd.concat([tab,df[sample]], axis = 1)
        i = 1
    del df

    # Re-format data table and merge with annotation data table
    tab['0'] = tab['0'].str.replace(' \|','')
    tab['0'] = tab['0'].str.strip()
    df = tab.transpose()
    df.reset_index(inplace=True)
    df.columns = df.iloc[0]
    df.drop(df.index[0], inplace = True)
    df.rename(columns={"0": "ID"}, inplace = True)

    # Calculated mapping rates
    tab = pd.merge(df, pheno, on=['ID'])
    tab['TOTAL_READS'] = pd.to_numeric(tab['Number of input reads'])
    tab['UNMAPPED_READS'] = pd.to_numeric(tab['Number of reads unmapped: too many mismatches']) + \
                            pd.to_numeric(tab['Number of reads unmapped: too short']) + \
                            pd.to_numeric(tab['Number of reads unmapped: other']) + \
                            pd.to_numeric(tab['Number of chimeric reads'])
    tab['MAPPED_READS'] = tab['TOTAL_READS'] - tab['UNMAPPED_READS']
    tab['MAPPED_IN_PERC'] = tab['MAPPED_READS'] / tab['TOTAL_READS']
    tab['UNMAPPED_IN_PERC'] = tab['UNMAPPED_READS'] / tab['TOTAL_READS']
    df = tab[['ID','GROUP','TOTAL_READS','MAPPED_READS','MAPPED_IN_PERC','UNMAPPED_READS','UNMAPPED_IN_PERC']]

    # Write table to output file
    df.to_csv(outfile, sep = '\t', index = False)


# --------------------------------------------------------------------
# Determine the git repository commit hash to create git url with hash
# e.g. https://github.roche.com/schmucr1/bksnake/commit/30a27ee
def get_git_revision_short_hash() -> str:
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()



# ------------------------------------------------------------------------------
# Define helper functions
# For counting number of lines in (fastq) file
def count_gzip_lines(filename):
    if not os.path.isfile(filename):
        raise Exception('count_gzip_lines: file is not present:'+filename)
    elif os.path.getsize(filename) == 0:
        raise Exception('count_gzip_lines: file is empty, i.e. getsize is zero:'+filename)
    else:
        with gzip.open(filename, 'rb') as f:
            for i, l in enumerate(f):
                pass
            #    print("File {1} contain {0} lines".format(i + 1, filename))
        return i+1
    

# ------------------------------------------------------------------------------
# Check whether the input sample IDs are unique or not    
def check_unique_sampleids(ids):
    if len(ids)!=len(set(ids)):
        seen = {}
        dups = []
        for id in ids:
            if id not in seen:
                seen[id]=1
            else:
                if seen[id] == 1:
                    dups.append(id)
                seen[id]+=1
        print('Duplicated ids:' + ', '.join(dups), file=sys.stderr)
        print('Duplicated sample ids detected. Program quits.', file=sys.stderr)
    else:
        pass

# ------------------------------------------------------------------------------
# Read sample annotations from a flat file
def read_biokit_sample_anno(filename):
    """read biokit sample annotation file, ensure that file
    paths are absolute paths, and returns a DataFrame"""

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


# ------------------------------------------------------------------------------
# Convert sample annotations to PhenoData format
def biokit_sample_anno_to_phenoData(anno):
    """convert biokit sample annotation DataFrame into the phenoData DataFrame
    in the same format as the biokit output pipeline.

    Args:
        anno: DataFrame returned by read_biokit_sample_anno

    Returns:
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


# ------------------------------------------------------------------------------
# Translate input biokit sample annotation file to output phenoData.meta
# used, for example, in workflow/Snakefile
def translate_biokit_to_phenoData_meta(infile, outfile):
    """Translate input biokit sample annotation file to output phenoData.meta
 
    Args:
        infile: input sample annotation file
        output: output phenoData.meta file
    """
    anno = read_biokit_sample_anno(infile)
    outdf = biokit_sample_anno_to_phenoData(anno)
    outdf.to_csv(outfile, sep='\t', index=False)