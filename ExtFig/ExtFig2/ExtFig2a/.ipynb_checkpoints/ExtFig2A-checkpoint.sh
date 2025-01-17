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




Rscript ExtFig2A.R 