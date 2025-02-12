#/bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y

# Activate deeptools
conda activate motifbreakR

Rscript motifbreakR_HAR.R