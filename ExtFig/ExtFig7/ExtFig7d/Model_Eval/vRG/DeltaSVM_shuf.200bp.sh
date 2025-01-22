#bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y


perl /shen/shenlabstore3/ijones1/GKM_explain_test/MG/Delta/deltasvm_script/deltasvm.pl /shen/shenlabstore3/ijones1/GKM_explain_test/SCZ_snps/Scz_Ref_Shuffled.SNPs.200bp.5.1.24.V2.fa /shen/shenlabstore3/ijones1/GKM_explain_test/SCZ_snps/Scz_Alt_Shuffled.SNPs.200bp.5.1.24.V2.fa ../Full_training_update/KMER.11.vRG.Full.model.8.5.24.model.txt  Delta_Shuf_Scz.200bp.Filt.allATAC.txt