#bin/bash
#$ -l h_rt=48:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y
#$ -t 1-4
#$ -pe smp 8


INDEX=fasta.txt
fasta=$(awk 'FNR == i {print}' i=${SGE_TASK_ID} $INDEX)

INDEX2=output.txt
output=$(awk 'FNR == i {print}' i=${SGE_TASK_ID} $INDEX2)

model=/shen/shenlabstore3/ijones1/GKM_explain_test/oRG/Full_training_update/oRG.Full.model.8.5.24.model.txt


/shen/shenlabstore3/ijones1/dependencies/GKM_Explain/rlsgkm/bin/gkmpredict $fasta $model ISM.${output}

