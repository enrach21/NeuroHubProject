#/bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y

# Author: Ian Jones
# Date: 7.12.23
# Determine methylation levels within 2kb

source ~/.bash_profile
conda activate deeptools

ATAC=/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/PEAKS/Merged.4.11.23.overlap.optimal_peak.final.bed

Cell=MG

awk '$9=$2' ${Cell}_CG.report.10.cov.bed | awk 'OFS="\t" {print $1,$2-1,$9,$3,$4,$5}' > temp_${Cell}_CG.report.10.cov.bed 

bedtools intersect -a $ATAC -b temp_${Cell}_CG.report.10.cov.bed -wa -wb > ATAC_${Cell}_CpG.bismark.10.cov.bed


bedtools intersect -a $ATAC -b temp_${Cell}_CG.report.10.cov.bed -v > ATAC_${Cell}_no_CpG.bed

rm temp_${Cell}_CG.report.10.cov.bed 

awk '{print $0"\t1"}' ATAC_${Cell}_CpG.bismark.10.cov.bed > ATAC2_${Cell}_CpG.bismark.10.cov.bed

rm ATAC_${Cell}_CpG.bismark.10.cov.bed

bedtools merge -i ATAC2_${Cell}_CpG.bismark.10.cov.bed -c 8,9,10 -o sum > ATAC_binned_${Cell}_CpG.bismark.10.cov.bed

rm ATAC2_${Cell}_CpG.bismark.10.cov.bed

# add the ratio of methylation 
awk '{$7 = $4 / ($4 + $5)}1' ATAC_binned_${Cell}_CpG.bismark.10.cov.bed > ATAC_binned2_${Cell}_CpG.bismark.10.cov.bed

rm ATAC_binned_${Cell}_CpG.bismark.10.cov.bed
