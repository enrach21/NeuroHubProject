#/bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y

date
hostname

# activate enviroment:
source ~/.bashrc
conda activate homer

# Names
IN_DIR=/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/Intervene_Both/sets
BED=0001_MG.bed 
OUTPUT=/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/Homer
DIR=MG_6_27_24

mkdir $DIR



findMotifsGenome.pl ${IN_DIR}/${BED} hg38 ${OUTPUT}/${DIR} -size given
