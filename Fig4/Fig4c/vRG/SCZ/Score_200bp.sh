#bin/bash
#$ -l h_rt=48:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y
#$ -t 1-4
#$ -pe smp 8


INDEX=/shen/shenlabstore3/ijones1/GKM_explain_test/vRG/SCZ_score_update/fasta.200bp.txt
fasta=$(awk 'FNR == i {print}' i=${SGE_TASK_ID} $INDEX)

INDEX2=/shen/shenlabstore3/ijones1/GKM_explain_test/vRG/SCZ_score_update/output.200bp.txt
output=$(awk 'FNR == i {print}' i=${SGE_TASK_ID} $INDEX2)

model=/shen/shenlabstore3/ijones1/GKM_explain_test/vRG/Full_training_update/vRG.Full.model.8.5.24.model.txt


/shen/shenlabstore3/ijones1/dependencies/GKM_Explain/rlsgkm/bin/gkmexplain $fasta $model $output

