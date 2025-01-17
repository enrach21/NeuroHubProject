#/bin/bash
#$ -l h_rt=4:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=200G
#$ -j y
#$ -r y

awk '$9=$2' MG_CG.report.10.cov.bed  | awk '{print $1,$2-1,$9,$8}'  | sort -k1,1 -k2,2n | awk -v OFS='\t' '{ $1=$1; print}' > MG_CG.report.10.cov.bedgraph

grep -v chrMT MG_CG.report.10.cov.bedgraph > MG_CG.report.10.cov.filt.bedgraph 

/shen/shenlabstore3/ijones1/dependencies/bedGraphToBigWig MG_CG.report.10.cov.filt.bedgraph  /shen/shenlabstore3/shared/reference_genome/hg38/hg38.chrom.sizes MG_CG.report.10.cov.bw

