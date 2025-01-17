#bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y


perl /shen/shenlabstore3/ijones1/GKM_explain_test/MG/Delta/deltasvm_script/deltasvm.pl /shen/shenlabstore3/ijones1/GKM_explain_test/MG_AD_Snps/MG.AD.Filt.allATAC.SNPs.200bp.V2.fa /shen/shenlabstore3/ijones1/GKM_explain_test/MG_AD_Snps/MG.AD.Filt.allATAC.ALT.SNPs.200bp.V2.fa ../Full_training_update/KMER.11.MG.Full.model.8.5.24.model.txt Delta_AD.Filt.MG.ATAC.200bp.txt