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
CIBER=/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/CiberSort/Results/CIBERSORTx_Job12_Results.csv
# Chosen Marker genes
META=/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/RNA-seq/CiberSort/NH.RNA.Label.csv


Rscript ExtFig1E.R $CIBER $META