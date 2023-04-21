#!/bin/bash

function usage()
{
  echo ""
  echo "1 input arguments required:"
  echo "1. input config.yaml file"
  echo ""
  echo "2 optional arguments"
  echo "1. target rule for local run"
  echo "2. nr of cores for local run"
  echo "Contact roland.schmucki@roche.com / tel 71330"
  echo ""
  echo ""
}

get_abs_filename() {
  # $1 : relative filename
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

if [[ $# -lt 1 ]];then usage && return; fi

SNAKE_FILE=workflow/Snakefile_public
CONFIG_FILE=$(get_abs_filename $1)
GRAPH_FILE=bk-public_rulegraph.pdf
LSF_PROFILE=/apps/rocs/etc/apps/snakemake/lsf/v1.4_memfix

if [ ! -e $CONFIG_FILE ] ; then usage && return; fi

user=`whoami`
homeDir=/home/${user}

# check if arvados variables are set
ARVADOS_CONF="$HOME/.credentials/arkau.conf"
ARVADOS_API_HOST=${ARVADOS_API_HOST:-$(grep "ARVADOS_API_HOST=" "$ARVADOS_CONF" | \
             sed -e "s/^export //" | cut -f2 -d=)}
ARVADOS_API_TOKEN=${ARVADOS_API_TOKEN:-$(grep "ARVADOS_API_TOKEN=" "$ARVADOS_CONF"  | \
             sed -e "s/^export //" | cut -f2 -d=)}

if [[ -z "$ARVADOS_API_HOST" || -z "$ARVADOS_API_TOKEN" ]]; then
  echo *Warning*: ARVADOS_API_TOKEN or ARVADOS_API_HOST variables are not set.
  echo Try setting them via ~/.bashrc or $ARVADOS_CONF file.
fi
export SINGULARITYENV_ARVADOS_API_HOST="$ARVADOS_API_HOST"
export SINGULARITYENV_ARVADOS_API_TOKEN="$ARVADOS_API_TOKEN"

echo "Load modules..."
ml purge && ml snakemake

if [[ $# -gt 1 ]]; then
  TARGET_RULE=$2
  if [[ $# -gt 2 ]]; then
    CORES=$3
  else
    CORES=1
  fi
  echo "Launch pipeline..."
  snakemake --snakefile $SNAKE_FILE \
       --configfile $CONFIG_FILE \
       --cores $CORES \
       --notemp \
       --latency-wait 6 \
       --rerun-incomplete \
       --keep-going \
       --verbose \
       --use-singularity \
       --singularity-prefix /projects/site/pred/ngs/Singularity/images \
       --singularity-args "--contain --cleanenv \
                           --env AWS_CONFIG_FILE=/projects/site/pred/ngs/pipelines/bksnake/.aws/config \
                           --env AWS_SHARED_CREDENTIALS_FILE=/projects/site/pred/ngs/pipelines/bksnake/.aws/credentials \
                           --bind /tmp \
                           --bind /projects/site/pred/ngs \
                           --bind '${homeDir}'/.credentials:/credentials \
                           --bind /projects/site/pred/ngs/pipelines/bksnake/.aws \
                           --bind /projects/site/pred/ngs/pipelines/bksnake/.credentials \
                           --bind /projects/site/pred/ngs/pipelines/bksnake/.credentials/ribiosAnnotation-secrets.json:/root/.credentials/ribiosAnnotation-secrets.json \
                           --bind /apps/rocs/2020.08/prefix/etc/ssl/certs \
                           --env REQUESTS_CA_BUNDLE=/apps/rocs/2020.08/prefix/etc/ssl/certs/ca-certificates.crt \
                           --env SSL_CERT_FILE=/apps/rocs/2020.08/prefix/etc/ssl/certs/ca-certificates.crt" ${TARGET_RULE}
else
  echo "Create graph..."
  snakemake --rulegraph all --snakefile $SNAKE_FILE --configfile $CONFIG_FILE | dot -Tpdf > $GRAPH_FILE

  echo "Launch pipeline..."
  snakemake --snakefile $SNAKE_FILE  \
       --configfile $CONFIG_FILE \
       --jobs 200 \
       --profile $LSF_PROFILE \
       --latency-wait 6 \
       --rerun-incomplete \
       --keep-going \
       --verbose \
       --use-singularity \
       --singularity-prefix /projects/site/pred/ngs/Singularity/images \
       --singularity-args "--contain --cleanenv \
                           --env AWS_CONFIG_FILE=/projects/site/pred/ngs/pipelines/bksnake/.aws/config \
                           --env AWS_SHARED_CREDENTIALS_FILE=/projects/site/pred/ngs/pipelines/bksnake/.aws/credentials \
                           --bind /tmp \
                           --bind /projects/site/pred/ngs \
                           --bind '${homeDir}'/.credentials:/credentials \
                           --bind /projects/site/pred/ngs/pipelines/bksnake/.aws \
                           --bind /projects/site/pred/ngs/pipelines/bksnake/.credentials \
                           --bind /projects/site/pred/ngs/pipelines/bksnake/.credentials/ribiosAnnotation-secrets.json:/root/.credentials/ribiosAnnotation-secrets.json \
                           --bind /apps/rocs/2020.08/prefix/etc/ssl/certs \
                           --env REQUESTS_CA_BUNDLE=/apps/rocs/2020.08/prefix/etc/ssl/certs/ca-certificates.crt \
                           --env SSL_CERT_FILE=/apps/rocs/2020.08/prefix/etc/ssl/certs/ca-certificates.crt"
fi
