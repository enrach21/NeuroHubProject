#bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y


# Activate Enviroment
source ~/.bashrc
conda activate r_analysis

python /shen/shenlabstore3/ijones1/dependencies/MEME/meme-5.5.5/scripts/fasta-dinucleotide-shuffle.py -f Scz.Filt.allATAC.SNPs.200bp.5.1.24.fa -c 3 > Scz.Shuffled.SNPs.200bp.5.1.24.fa