#! python

import sys
import pandas as pd
import os.path
import numpy as np



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

# original code
#    id_group = anno.iloc[:,0].astype(str) + "_" + anno.iloc[:,1].astype(str)
#    keepcols = [True] * len(anno.columns)
#    keepcols[2:4] = [False]*2 # FASTQ1 and FASTQ2
#    rest = anno.iloc[:, keepcols].copy()
#    res = pd.concat([id_group, rest], axis=1)
#    res.rename(columns={res.columns[0]: "ID_GROUP",
#        res.columns[1]: "ID",
#        res.columns[2]: "GROUP"
#        }, inplace=True)

# modified code:
# Do NOT concatenate ID and GROUP - keep only ID (for test)
#    print("biokit_sample_anno_to_phenoData:")
#    print(anno)
    anno['GROUP'] = anno['GROUP'].astype(str)
    anno['GROUP'] = anno['GROUP'].replace(to_replace=r'^NA$', value='noGroup', regex=True)

    id_group = anno.iloc[:,0].astype(str) + "_" + anno.iloc[:,1].astype(str)
    keepcols = [True] * len(anno.columns)
    keepcols[2:4] = [False]*2 # FASTQ1 and FASTQ2
    rest = anno.iloc[:, keepcols].copy()
    res = pd.concat([rest, id_group], axis=1)
    res.rename(columns={'#ID': 'ID', 0:'ID_GROUP'}, inplace=True)

    return(res)


def translate_biokit_to_phenoData_meta(infile, outfile):
    """Translate input biokit sample annotation file to output phenoData.meta
 
    Args:
        infile: input sample annotation file
        output: output phenoData.meta file
    """
    anno = read_biokit_sample_anno(infile)
    outdf = biokit_sample_anno_to_phenoData(anno)
    outdf.to_csv(outfile, sep='\t', index=False)



def fastqc_filename(filename):
    """return fastqc file name of a FASTQ file"""
    fdir = os.path.dirname(filename)
    bn = os.path.basename(filename)
    suffind = bn.find('.f')
    if suffind == -1:
        raise(Exception("the file " + bn + " does look like FASTQ file"))
    outbn = bn[0:suffind]+'_fastqc.zip'

    return(outbn)



if __name__ == '__main__':
    res = read_biokit_sample_anno("../data/sample_annotation.txt")
    sampleid = list(res.iloc[:, 0])
    fastqs = list(res.iloc[:, 2])
    for f in fastqs:
        print(f)
        print("-->")
    phenoData = biokit_sample_anno_to_phenoData(res)
    print(phenoData)
