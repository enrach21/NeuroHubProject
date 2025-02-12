#bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y


perl /shen/shenlabstore3/ijones1/GKM_explain_test/MG/Delta/deltasvm_script/deltasvm.pl /shen/shenlabstore3/ijones1/GKM_explain_test/HARs/Ref_Shuffled.HAR.allATAC.variants.7.22.24.V2.fa /shen/shenlabstore3/ijones1/GKM_explain_test/HARs/Alt_Shuffled.HAR.allATAC.variants.7.22.24.V2.fa  ../Full_training_update/KMER.11.oRG.Full.model.8.5.24.model.txt Shuf_Delta_HAR_Score_7_22_24.txt