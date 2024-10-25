#bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y


perl /shen/shenlabstore3/ijones1/GKM_explain_test/MG/Delta/deltasvm_script/deltasvm.pl /shen/shenlabstore3/ijones1/GKM_explain_test/SCZ_snps/Scz.Filt.allATAC.SNPs.5.1.24.V2.fa /shen/shenlabstore3/ijones1/GKM_explain_test/SCZ_snps/Scz.Filt.allATAC.ALT.SNPs.5.1.24.V2.fa ../Full_training_update/KMER.11.MG.Full.model.8.5.24.model.txt  Delta_Scz.Filt.allATAC.txt