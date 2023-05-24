"""
Wrapper script for the "bulk RNASeq biokit" Snakemake workflow
"""
import os
import sys
import yaml
import logging
import argparse
import subprocess
from pathlib import Path


"""
Declare global variables
"""
GRAPHFILE = 'rulegraph.pdf' # name of the output graph file
GRAPHFILE_PNG = 'rulegraph.png'
CONFIGFILE = 'config.yaml'  # name of the output config file
CONFIG = None # global variable for input config dictionary


"""
Verify whether the provided command/path for Snakemake is working
"""
def verify_snakemake(path):
    cmd = path+' --version'
    res = subprocess.run(cmd, shell=True, check=False, env=os.environ.copy(), stdout=subprocess.DEVNULL)
    if res.returncode != 0:
        raise Exception ('No executable found for Snakemake. Please specify the parameter \'snakemake path\'')


"""
Read common configuration file
Normally, this is file 'config/common.yaml'
"""
def read_config(file):
    yaml.preserve_quotes = True

    if not os.path.isfile(file):
        raise Exception('Common configuration file does not exist. Abort! '+file)
 
    res = yaml.safe_load(Path(file).read_text())
    return res


"""
Create the config file for the Snakemake run

Use the input config file (default: config/common.yaml) and 
add parameters from the user input. CONFIG is a global variable.
Then, output the as a new file 'config.yaml' into the
output directory

Note that these 3 variables are not defined in the default config:
CONFIG['species']
CONFIG['species_name']
CONFIG['organism']
"""
def create_run_config(args):

    # update config dictionary
    if args.outdir:
        CONFIG['outdir'] = args.outdir

    if args.genome_dir:
        CONFIG['genome_dir'] = args.genome_dir
        
    if args.species:
        CONFIG['species'] = CONFIG['genomes'][args.species]['species']
        CONFIG['species_name'] = CONFIG['genomes'][args.species]['species_name']
        CONFIG['organism'] = CONFIG['genomes'][args.species]['organism']
    
    if args.keep_fastq:
        CONFIG['keep_fastq_files'] = True
    else:
        CONFIG['keep_fastq_files'] = False

    if args.keep_bam:
        CONFIG['keep_bam_files'] = True
    else:
        CONFIG['keep_bam_files'] = False

    if args.not_generate_cram:
        CONFIG['generate_cram_files'] = False
    else:
        CONFIG['generate_cram_files'] = True

    if args.generate_bw:
        CONFIG['generate_bw_files'] = True
    else:
        CONFIG['generate_bw_files'] = False

    if args.generate_unmapped:
        CONFIG['generate_unmapped'] = True
    else:
        CONFIG['generate_unmapped'] = False
        
    if args.cutadapt:
        CONFIG['cutadapt']['run'] = True
        if args.cutadapt_params:
            CONFIG['cutadapt']['parameters'] = args.cutadapt_params
    else:
        CONFIG['cutadapt']['run'] = False
        
    if args.snakemake_path:
        CONFIG['snakemake']['path'] = args.snakemake_path

    if args.snakemake_parameters:
        CONFIG['snakemake']['parameters'] = args.snakemake_parameters
        
    if args.singularity_prefix:
        CONFIG['singularity']['prefix'] = args.singularity_prefix

    if args.metadata_file:
        CONFIG['metadata']['file'] = args.metadata_file
    if args.metadata_group_name:
        CONFIG['metadata']['group_name'] = args.metadata_group_name

    # check output config file
    configfile = os.path.join(args.outdir, CONFIGFILE)
#    if os.path.isfile(configfile):
#        sys.stderr.write('Configuration file exists already: '+configfile)

    # write dictionary tp new output yaml file
    f = open(configfile, 'w')
    yaml.dump(CONFIG, f)
    f.close()
   
    
"""
Create Snakemake Graph file
"""
def graph(args):
    snakefile = os.path.abspath(args.snakefile)
    outfile = os.path.join(args.outdir, GRAPHFILE)
    outfile2 = os.path.join(args.outdir, GRAPHFILE_PNG)
    configfile = os.path.join(args.outdir, CONFIGFILE)
                        
    cmd = CONFIG['snakemake']['path']+' '\
          +' --rulegraph all'\
          +' --snakefile '+snakefile\
          +' --configfile '+configfile\
          +' | dot -Tpdf > '+outfile\
          +' && pdftoppm -png '+outfile+' > '+outfile2

    print('\n=======================================\n'+
          'Command to re-run the pipeline-graph:\n'+cmd+
          '\n=======================================\n\n')
    res = subprocess.run(cmd, shell=True, check=True, env=os.environ.copy())
    print('exit status code:', res.returncode )
    return res
    
    
"""
Run Snakemake workflow
  --notemp, --nt   Ignore temp() declarations. This is useful when running only a part of the workflow, 
                   since temp() would lead to deletion of probably needed files by other parts of the workflow.
    
  --jobs [N], -j [N]    Use at most N CPU cluster/cloud jobs in parallel. For local execution this 
                        is an alias for --cores. (default: None)
"""
def workflow(args):
    snakefile = os.path.abspath(args.snakefile)
    configfile = os.path.join(args.outdir, CONFIGFILE)
    
    homeDir = str(Path.home())

    cmd = CONFIG['snakemake']['path']+' '\
          +' --snakefile '+snakefile\
          +' --configfile '+configfile\
          +' --latency-wait 6'\
          +' --rerun-incomplete'\
          +' --keep-going'\
          +' --verbose'\
          +' '+CONFIG['snakemake']['parameters']\
          +' --use-singularity'\
          +' --singularity-prefix '+CONFIG['singularity']['prefix']\
          +' --singularity-args "--contain --cleanenv'\
          +' --bind '+CONFIG['genome_dir']\
          +' --bind /tmp"'
    
    if args.jobs:
        cmd = cmd+' --jobs '+str(args.jobs)+' --profile '+args.cluster_profile
        print('Submit workflow to the cluster')
    else:
        cmd = cmd+' --cores '+str(args.cores)
        print('Run workflow locally')
    
    if args.target:
        cmd = cmd+' \''+args.target+'\''

    print('\n=======================================\n'+
          'Command to re-run the pipeline:\n'+cmd+
          '\n=======================================\n\n')
    res = subprocess.run(cmd, shell=True, check=True, env=os.environ.copy())
    print('exit status code:', res.returncode )
    return res



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Wrapper script for the Snakemake workflow of the biokit bulk RNASeq pipeline aka bksnake.'
    )

    # Required arguments
    parser.add_argument('--outdir', '-o', dest='outdir', 
                        help='Path to output directory. Name in config: \'outdir\'', 
                        required=True, default=None)

    # Optional arguments
    # Parameters that are also present in the Snakemake config file
    # All default values from that config file (template or user given)
    parser.add_argument('--metadata-file',
                        help='Path to input tab-delimited text file containing metadata information. Name in config: \'metadata: file:\'', 
                        required=False, default=None)
    parser.add_argument('--metadata-group-name', dest='metadata_group_name',
                        help='Metadata column name that should be used as condition category for the qc analysis (e.g. pca plot) \
                        in case of combining multiple columns, then use semi-colon as separator, e.g. \'Tissue;Treatment ID\'. Name in config: \'metadata: group_name\'',
                        required=False, default=None)
    parser.add_argument('--genome-dir',
                        help='Path to the genome root directory. Must contain sub-directory for each species genomes. Name in config: \'genome_dir\'', 
                        required=False, default=None)
    parser.add_argument('--singularity-prefix',
                        help='Path to the location where Singularity images should be stored for re-use later. Name in config: \'singularity: prefix\'', 
                        required=False, default=None)
    parser.add_argument('--keep-fastq',
                        help='Keep a copy of the input fastq files in the output directory. Name in config: \'keep_fastq_files\'', 
                        required=False, action='store_true')
    parser.add_argument('--keep-bam',
                        help='Keep bam files in the output directory. Name in config: \'keep_bam_files\'', 
                        required=False, action='store_true')
    parser.add_argument('--not-generate-cram',
                        help='Do generate cram files. Name in config: \'generate_cram_files\'', 
                        required=False, action='store_true')
    parser.add_argument('--generate-bw',
                        help='Generate bigwig, aligned reads coverage files. Name in config: \'generate_bw_files\'', 
                        required=False, action='store_true')
    parser.add_argument('--generate-unmapped',
                        help='Generate files with unmapped reads. Name in config: \'generate_unmapped\'', 
                        required=False, action='store_true')
    parser.add_argument('--cutadapt',
                        help='Run Cudadapt read trimming prior to read mapping. Name in config: \'cutadapt: run\'', 
                        required=False, action='store_true')
    parser.add_argument('--cutadapt-parameters', dest='cutadapt_params',
                        help='Additional Cutadapt parameters, e.g. --cutadapt-parameters=\'--minimum-length 10\'. Name in config: \'cutadapt: parameters\'', 
                        required=False, default=None)
    parser.add_argument('--snakemake-path',
                        help='Path to Snakemake. Name in config: \'snakemake: path\'',
                        required=False, default=None)
    parser.add_argument('--snakemake-parameters',
                        help='Additional Snakemake parameters, e.g. --snakemake-parameters=\'--dry-run --notemp\'. Name in config: \'snakemake: parameters\'',
                        required=False, default=None)



    # Parameters not present in Snakemake config file 
    parser.add_argument('--species', '-s', dest='species',
                        help='Reference genome species for mapping, e.g. hg38, mm39, mfa5, rn7, ss11, oc2. Not in config.',
                        required=False, default=None)
    parser.add_argument('--config', '-f',
                        help='Path to input yaml config file for Snakemake. All parameters of the config file are \
                        overwritten if they are specified by optional arguments to this wrapper script',
                        required=False, default='config/config.yaml')
    parser.add_argument('--snakefile',  
                        help='Path to input Snakemake file.', 
                        required=False, default='workflow/Snakefile')
    parser.add_argument('--target', '-t', 
                        help='Target rule for Snakemake workflow, e.g. all', 
                        required=False, default='all')
    parser.add_argument('--cluster-profile', '-p', 
                        help='Path to cluster profile, only used with --jobs', 
                        required=False, default='config/lsf')
    parser.add_argument('--no-dag',
                        help='Do not create \"directed acyclic graph\" of the workflow.', 
                        required=False, action='store_true')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--cores', '-c', help='Number of cores to use for local run, required for local run, e.g. 8')
    group.add_argument('--jobs', '-j', help='Number of jobs for running on the cluster, required for run on cluster, e.g. 100')
    
    args = parser.parse_args()

    # Check snakemake file
    if not os.path.isfile(args.snakefile):
        raise Exception('Snakemake file does not exist. Abort! '+args.snakefile)
    
    # Check input or (default:template) config file
    if not os.path.isfile(args.config):
        raise Exception('Configuration (template) file does not exist. Abort! '+args.config)
   
    # Create output directory
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    
    # Read input config file and return as global variable CONFIG
    CONFIG = read_config(args.config)

    # Get species: 
    # from command line -> now required.
    # from config file
    # from metadata file
    
    # Create config.yaml file
    create_run_config(args)
 
    # Verify Snakemake path is executable
    verify_snakemake(CONFIG['snakemake']['path'])

    # Create graph file
    if args.no_dag == False:
        res = graph(args)
    
    # Launch workflow
    res = workflow(args)
