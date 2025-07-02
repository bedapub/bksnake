#!/bin/bash

log=$1
input_log2=$2
input_pheno=$3
input_collapsed=$4
output_cls=$5
output_tab=$6
output_txt=$7
output_thr=$8
output_pca=$9
output_bioqc=${10}

R -e "df=ribiosIO::readTable('${input_pheno}'); ribiosIO::write_cls(factor(df\$GROUP),'${output_cls}')" > $log

(plotPCA.Rscript -infile ${input_log2} -outfile ${output_pca} -outTable ${output_tab} -cls ${output_cls}) >> $log

(bioqc.Rscript -featuretype GeneSymbol -infile ${input_collapsed} -outfile ${output_txt}) >> $log

(biosHeatmap.Rscript -infile ${output_txt} -outfile ${output_bioqc} -scale none -colors blackred \
  -naColor lightgray -symbreaks auto -dendrogram both -dist euclidean -hclust ward.D2 \
  -xlab -ylab -cexRow -cexCol -colorKeyTitle -width -height -margins -zlimLo -zlimHi -main 'BioQC') >> $log

(bioqc.Rscript -threshold 2.0 -featuretype GeneSymbol -infile ${input_collapsed} -outfile ${output_thr}) >> $log
