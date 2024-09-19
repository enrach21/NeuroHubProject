#bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y


# Activate Enviroment
source ~/.bashrc


python /shen/shenlabstore3/ijones1/dependencies/MEME/meme-5.5.5/scripts/fasta-dinucleotide-shuffle.py -f MG.AD.Filt.allATAC.SNPs.200bp.fa -c 3 > Shuffled.MG.AD.Filt.allATAC.SNPs.200bp.fa