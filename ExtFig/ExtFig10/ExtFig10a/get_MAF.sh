#/bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -m ea
#$ -M Ian.Jones3@ucsf.edu

# Activate deeptools
conda activate MAFsinRegion

mafsInRegion oRG.ATAC.HARs.bed -outDir  /shen/shenlabstore3/ijones1/dependencies/vcf2maf/HAR_MAF_V2 /shen/shenlabstore3/ijones1/dependencies/vcf2maf/MAF/*