#bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y

# Author: Ian Jones, Ian.Jones3@ucsf.edu

# Dependencies please make sure that you have bowtie2, bedtools and samtools installed.
source ~/.bash_profile
conda activate r_analysis



#### Get barplot for PR0.9 ####
# Location of normalized expression
RNA=/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/Marker_expression/AVG_TMM_RPKM_NoLog_exp_geneID_9.20.23.csv
# Chosen Marker genes
MARKER=/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/NH.Select.Markers.V2.csv


Rscript Fig1C.R $RNA $MARKER 