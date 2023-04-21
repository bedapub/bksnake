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
CONFIGFILE = 'config.yaml' # name of the output config file
COMMON_CONFIG = None # global variable for input common config file


"""
Read common configuration file
"""
def read_common_config(file):
    yaml.preserve_quotes = True

    if not os.path.isfile(file):
        raise Exception('Common configuration file does not exist. Abort! '+file)
 
    res = yaml.safe_load(Path(file).read_text())
    return res


"""
Create the config file for the Snakemake run

Use template and add parameters from user to the output config.yaml

template: config/config_template.yaml
"""
def create_run_config(args):

    yaml.preserve_quotes = True

    # read input config file (template)
    conf = yaml.safe_load(Path(args.config).read_text())

    # update config dictionary
    if args.outdir:
        conf['results'] = args.outdir
        
    if args.species:
        conf['species'] = COMMON_CONFIG['genomes'][args.species]['species']
        conf['species_name'] = COMMON_CONFIG['genomes'][args.species]['species_name']
        conf['organism'] = COMMON_CONFIG['genomes'][args.species]['organism']
        conf['genome_dir'] = COMMON_CONFIG['genomes'][args.species]['genome_dir']
    
    if args.keep_fastq:
        conf['keep_fastq_files'] = True
    else:
        conf['keep_fastq_files'] = False

    if args.keep_bam:
        conf['keep_bam_files'] = True
    else:
        conf['keep_bam_files'] = False

    if args.not_generate_cram:
        conf['generate_cram_files'] = False
    else:
        conf['generate_cram_files'] = True

    if args.generate_bw:
        conf['generate_bw_files'] = True
    else:
        conf['generate_bw_files'] = False

    if args.generate_unmapped:
        conf['generate_unmapped'] = True
    else:
        conf['generate_unmapped'] = False
        
    if args.cutadapt:
        conf['cutadapt']['run'] = True
        if args.cutadapt_params:
            conf['cutadapt']['parameters'] = args.cutadapt_params
    else:
        conf['cutadapt']['run'] = False

    # check output config file
    configfile = os.path.join(args.outdir, CONFIGFILE)
    if os.path.isfile(configfile):
        #raise Exception('Configuration file exists already. Abort! '+configfile)
        sys.stderr.write('Configuration file exists already: '+configfile)

    # write dictionary tp new output yaml file
    f = open(configfile, 'w')
    yaml.dump(conf, f)
    f.close()
   
    
"""
Create Snakemake Graph file
"""
def graph(args):
    snakefile = os.path.abspath(args.snakefile)
    outfile = os.path.join(args.outdir, GRAPHFILE)
    configfile = os.path.join(args.outdir, CONFIGFILE)
                        
    cmd = COMMON_CONFIG['snakemake']+' '\
          +' --rulegraph all'\
          +' --snakefile '+snakefile\
          +' --configfile '+configfile\
          +' | dot -Tpdf > '+outfile
    
    print('\n=======================================\n'+
          'Command to re-run the pipeline-graph:\n'+cmd+
          '\n=======================================\n\n')
    res = subprocess.run(cmd, shell=True, check=True, env=os.environ.copy())
    print('exit status code:', res.returncode )
    return res
    
    
"""
Run Snakemake workflow
  --notemp, --nt   Ignore temp() declarations. This is useful when running only a part of the workflow, since temp() would lead to deletion of probably needed files by other parts of the workflow.
    
  --jobs [N], -j [N]    Use at most N CPU cluster/cloud jobs in parallel. For local execution this is an alias for --cores. (default: None)
"""
def workflow(args):
    snakefile = os.path.abspath(args.snakefile)
    configfile = os.path.join(args.outdir, CONFIGFILE)
    
    homeDir = str(Path.home())

    cmd = COMMON_CONFIG['snakemake']+' '\
          +' --snakefile '+snakefile\
          +' --configfile '+configfile\
          +' --latency-wait 6'\
          +' --rerun-incomplete'\
          +' --keep-going'\
          +' --verbose'\
          +' '+args.snake_params\
          +' --use-singularity'\
          +' --singularity-prefix '+COMMON_CONFIG['singularity_prefix']\
          +' --singularity-args "--contain --cleanenv'\
          +' --bind '+COMMON_CONFIG['genomes'][args.species]['genome_dir']\
          +' --bind /tmp"'
    
    if args.jobs:
        cmd = cmd+' --jobs '+str(args.jobs)+' --profile '+args.profile
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
                        help='Path to output directory', 
                        required=True, default=None)
    
    # Optional arguments
    parser.add_argument('--species', '-s', dest='species',
                        help='Reference genome species for mapping, e.g. hg38, mm39, mfa5, rn7, ss11, oc2', 
                        required=False, default=None)
    parser.add_argument('--config', '-f', 
                        help='Path to input yaml config file for Snakemake. All parameters of the config file are \
                        overwritten if they are specified by optional arguments to this wrapper script', 
                        required=False, default='config/config_template.yaml')
    parser.add_argument('--common-config', dest='common',
                        help='Path to common yaml config file for Snakemake.', 
                        required=False, default='config/common.yaml')
    parser.add_argument('--keep-fastq',
                        help='Keep a copy of the input fastq files in the output directory.', 
                        required=False, action='store_true')
    parser.add_argument('--keep-bam',
                        help='Keep bam files in the output directory.', 
                        required=False, action='store_true')
    parser.add_argument('--not-generate-cram',
                        help='Do generate cram files.', 
                        required=False, action='store_true')
    parser.add_argument('--generate-bw',
                        help='Generate bigwig, aligned reads coverage files.', 
                        required=False, action='store_true')
    parser.add_argument('--generate-unmapped',
                        help='Generate files with unmapped reads.', 
                        required=False, action='store_true')
    parser.add_argument('--cutadapt',
                        help='Run Cudadapt read trimming prior to read mapping.', 
                        required=False, action='store_true')
    parser.add_argument('--cutadapt-parameters', dest='cutadapt_params',
                        help='Additional Cutadapt parameters, e.g. --cutadapt-parameters=\'--minimum-length 10\'', 
                        required=False, default='')    
    parser.add_argument('--snakefile',  
                        help='Path to input Snakemake file.', 
                        required=False, default='workflow/Snakefile')
    parser.add_argument('--target', '-t', 
                        help='Target rule for Snakemake workflow, e.g. all', 
                        required=False, default='')
    parser.add_argument('--profile', '-p', 
                        help='Path to cluster profile, only used with --jobs', 
                        required=False, default='/apps/rocs/etc/apps/snakemake/lsf/v1.4_memfix')
    parser.add_argument('--no-dag',
                        help='Do not create \"directed acyclic graph\" of the workflow.', 
                        required=False, action='store_true')
    parser.add_argument('--snakemake-parameters', dest='snake_params',
                        help='Additional Snakemake parameters, e.g. --snakemake-parameters=\'--dry-run --notemp\'', 
                        required=False, default='')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--cores', '-c', help='Number of cores to use for local run, required for local run, e.g. 8')#, default=8)
    group.add_argument('--jobs', '-j', help='Number of jobs for running on the cluster, required for run on cluster, e.g. 100')#, default=100)
    
    args = parser.parse_args()

    # Check snakemake file
    if not os.path.isfile(args.snakefile):
        raise Exception('Snakemake file does not exist. Abort! '+args.snakefile)
    
    # Check input or (default:template) config file
    if not os.path.isfile(args.config):
        raise Exception('Configuration (template) file does not exist. Abort! '+args.config)

    # Read common config file (global variable COMMON_CONFIG)
    COMMON_CONFIG = read_common_config(args.common)
    
    # Create output directory
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    
    # Create config.yaml file
    create_run_config(args)
  
    # Create graph file
    if args.no_dag == False:
        res = graph(args)
    
    # Launch workflow
    res = workflow(args)