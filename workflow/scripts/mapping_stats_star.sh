#!/bin/bash

if [[ $# -lt 2 || $# -gt 2 ]]; then
  echo ""
  echo "2  input arguments required:"
  echo "1. input folder with STAR Log.final.out files"
  echo "2. samples file (2 columns required: 1st col sample name, 2nd col sample group)"
  echo ""
  echo "Contact roland.schmucki@roche.com / tel 71330"
  echo ""
  exit 1
fi

umask 2
set -e
if [[ $BIOKIT_VERBOSE == 1 ]]; then
  set -x
fi
indir=$1
sampleFile=$2

### Check input files
if [ ! -d $indir ]; then
  echo "Input folder $indir does not exist. Abort!"
  exit 1
fi
if [ ! -e $sampleFile ]; then
  echo "Input sample file $sampleFile does not exist. Abort!"
  exit 1
fi

### Mapping stats
maxsamples=`grep -v ^# $sampleFile | grep -vw ID_GROUP | wc -l | cut -d\  -f1`
declare -a samples=(`grep -v ^# $sampleFile | grep -vw ID_GROUP | cut -f2`) 
declare -a groups=(`grep -v ^# $sampleFile | grep -vw ID_GROUP | cut -f3`)

# Header line
printf "#Description" > $indir/star_Log_final.txt
for (( n = 0 ; n < $maxsamples ; n++ )) do
  i=${samples[$n]}
  printf "\t$i" >> $indir/star_Log_final.txt
done
printf "\n" >> $indir/star_Log_final.txt

# Data
s=0
cut -f1 $indir/${samples[0]}_Log.final.out > $s
for (( n = 0 ; n < $maxsamples ; n++ )) do
  i=${samples[$n]}
  cut -f2 $indir/${samples[$n]}_Log.final.out > $i.tmp
  s=`echo $s $i.tmp`
done
paste $s | sed 's/ |//g' | grep -vw READS >> $indir/star_Log_final.txt
rm $s

# Mapping Stats
for (( n = 0 ; n < $maxsamples ; n++ )) do
  i=${samples[$n]}
  f=$indir/${samples[$n]}_Log.final.out
  t=`grep "Number of input reads" $f | awk '{print $NF}'`
  m1=`grep "Uniquely mapped reads number" $f | awk '{print $NF}'`
  m2=`grep "Number of reads mapped to multiple loci" $f | awk '{print $NF}'`
  m=`echo $m1 + $m2 | bc`
  m3=`grep "Number of reads mapped to too many loci" $f | awk '{print $NF}'`
  m=`echo $m + $m3 | bc`
  m4=`grep "Number of chimeric reads" $f | awk '{print $NF}'`
  m=`echo $m + $m4 | bc`
  u=`echo $t - $m | bc`
  printf "%s\t%s\t%d\t%d\n" $i ${groups[$n]} $t $u >> tmp$$
done
printf "ID\tGROUP\tTOTAL_READS\tMAPPED_READS\tMAPPED_IN_PERC\tUNMAPPED_READS\tUNMAPPED_IN_PERC\n"
awk 'BEG{FS="\t"}{if ($3==0){r1=0.; r2=0.}else{r1=($3-$4)/$3;r2=$4/$3};printf("%s\t%s\t%d\t%d\t%f\t%d\t%f\n",$1,$2,$3,$3-$4,r1,$4,r2)}' tmp$$
rm tmp$$
exit 0
