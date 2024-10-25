#bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y

/shen/shenlabstore3/ijones1/dependencies/GKM_Explain/rlsgkm/bin/gkmpredict /shen/shenlabstore3/ijones1/GKM_explain_test/deltasvm_script/kmer.11.txt oRG.Full.model.8.5.24.model.txt KMER.11.oRG.Full.model.8.5.24.model.txt