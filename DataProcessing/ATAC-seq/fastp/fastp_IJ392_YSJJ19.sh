#/bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$ -cwd
#$ -j y
#$ -r y
#$ -m ea
#$ -M Ian.Jones3@ucsf.edu

# Example: qsub -v READ1=R1.fastq.gz  -v READ2=R2.fastq.gz -v SAMPLE=IJ001 ./fastp.sh

date
hostname

# activate enviroment:
source ~/.bashrc
conda activate base

# Define folders and load dependencies.
output=/shen/shenlabstore3/ijones1/NeuroHub/ATAC-seq/fastp

READ1=/shen/shenlabstore/sequencing_data/YSJJ19/IJ392_S105_L003_R1_001.fastq.gz 
READ2=/shen/shenlabstore/sequencing_data/YSJJ19/IJ392_S105_L003_R2_001.fastq.gz 
trim=100
PE=true
SAMPLE=IJ392_YSJJ19

# echo
echo $READ1
echo $trim
echo $output
echo $SAMPLE

# Perform analysis
if $PE
then
/shen/shenlabstore3/ijones1/dependencies/fastp/fastp -i $READ1  -I $READ2   -o ${output}/${SAMPLE}.fastp.R1.fq.gz -O ${output}/${SAMPLE}.fastp.R2.fq.gz -b ${trim} -B ${trim} -p -h "${output}/${SAMPLE}.fastp.html" -j "${output}/${SAMPLE}.fastp.json"
else
	fastp -i $READ1 -o ${output}/${SAMPLE}.fastp.R1.fq.gz -p -t ${trim} -h "${output}/${SAMPLE}.fastp.html" -j "${output}/${SAMPLE}.fastp.json"
fi
