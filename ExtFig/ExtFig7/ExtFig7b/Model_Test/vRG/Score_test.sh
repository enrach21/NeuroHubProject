#bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y
#$ -t 1-2
#$ -pe smp 8


for cell in 'MG' 'OPC' 'oRG' 'vRG'

do


    INDEX=/shen/shenlabstore3/ijones1/GKM_explain_test/${cell}/Test_Model_Update/fasta.txt
    fasta=$(awk 'FNR == i {print}' i=${SGE_TASK_ID} $INDEX)

    INDEX2=/shen/shenlabstore3/ijones1/GKM_explain_test/${cell}/Test_Model_Update/output.txt
    output=$(awk 'FNR == i {print}' i=${SGE_TASK_ID} $INDEX2)

    model=/shen/shenlabstore3/ijones1/GKM_explain_test/vRG/Full_training_update/vRG.Not_chr2.model.8.5.24.model.txt


    /shen/shenlabstore3/ijones1/dependencies/GKM_Explain/rlsgkm/bin/gkmpredict $fasta $model $output

done

