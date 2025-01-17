#bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y


# Activate Enviroment
source ~/.bashrc
conda activate GKM_SVM

Rscript /shen/shenlabstore3/ijones1/GKM_explain_test/Microglia/Neg_90k/Get_negative.R