#!/bin/bash

INFILE=$1
STRANDNESS=$2

if [[ $# -lt 2 || $# -gt 2 ]]; then
  echo ""
  echo "Outputs the strandness parameter for Picard, FeatureCounts, and Biokit tools."
  echo ""
  echo "2  input arguments required:"
  echo "1. output file from infer_experiment.py"
  echo "2. strandness mode parameter: auto, ff, or fr."
  echo ""
  echo "Contact roland.schmucki@roche.com"
  echo ""
  exit 1
fi

if [[ $STRANDNESS == "auto" ]]; then
    FEATURECOUNTS=$(awk '/^Fraction/ { UND=$NF *100; getline ; FF=$NF * 100 ; getline ; FR=$NF *100; diff=FF-FR; if ( UND > 50 || ( diff < 0 ? -diff : diff ) < 50 ) { print "0" } else { if ( FF > FR ) { print "1" } else { print "2" }  } }' $INFILE )
    if [[ $FEATURECOUNTS == "0" ]]; then
        DEFAULT="undetermined"
        PICARD="NONE"
        BIOKIT="0"
    elif [[ $FEATURECOUNTS == "1" ]]; then
        DEFAULT="stranded"
        PICARD="FIRST_READ_TRANSCRIPTION_STRAND"
        BIOKIT="1"
    elif [[ $FEATURECOUNTS == "2" ]]; then
        DEFAULT="stranded"
        PICARD="SECOND_READ_TRANSCRIPTION_STRAND"
        BIOKIT="2"
    else
        echo "Error in rule strandness: FEATURECOUNTS value is not 0, 1, nor 2. Abort!" && exit 1
    fi
else
    if [[ $STRANDNESS == "ff" ]]; then
        DEFAULT="stranded"
        FEATURECOUNTS="1"
        PICARD="FIRST_READ_TRANSCRIPTION_STRAND"
        BIOKIT="1"
    elif [[ $STRANDNESS == "fr" ]]; then
        DEFAULT="stranded"
        FEATURECOUNTS="2"
        PICARD="SECOND_READ_TRANSCRIPTION_STRAND"
        BIOKIT="2"
    else
        DEFAULT="undetermined"
        FEATURECOUNTS="0"
        PICARD="NONE"
        BIOKIT="0"
    fi
fi

echo "strandness_mode=$STRANDNESS"
echo "default=$DEFAULT"
echo "picard=$PICARD"
echo "featurecounts=$FEATURECOUNTS"
echo "biokit=$BIOKIT"