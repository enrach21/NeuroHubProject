#/bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y
#$ -m ea
#$ -M Ian.Jones3@ucsf.edu

# Date: 4/19/22
# Goal: Merge 1D peak calls into consensus anchor set

date
hostname

#### activate enviroment ####
source ~/.bashrc
conda activate deeptools

# Variable
output=/wynton/group/shen/NeuroHub.7.23.21/Final/PLAC/Anchor/Merged
CHROM=/shen/shenlabstore3/shared/reference_genome/hg38/hg38.chrom.sizes
oligo_anchor=${output}/oligo.shrt.30mil_peaks.broadPeak
microglia_anchor=${output}/MG.shrt.30mil_peaks.broadPeak 
oRG_anchor=${output}/oRG.shrt.30mil_peaks.broadPeak 
vRG_anchor=${output}/vRG.shrt.30mil_peaks.broadPeak



cat ${oligo_anchor} >> ${output}/Anchor.Merged.9.7.23.bed
cat ${microglia_anchor} >> ${output}/Anchor.Merged.9.7.23.bed
cat ${oRG_anchor} >> ${output}/Anchor.Merged.9.7.23.bed
cat ${vRG_anchor} >> ${output}/Anchor.Merged.9.7.23.bed


# Cut first 4 columns and resort
cut -f1-4 ${output}/Anchor.Merged.9.7.23.bed > ${output}/temp.bed 
sort -k1,1 -k2,2n ${output}/temp.bed > ${output}/Anchor.Merged.9.7.23.sorted.bed



# Merge any overlapping peaks
bedtools merge -i ${output}/Anchor.Merged.9.7.23.sorted.bed > ${output}/Anchor.Merged.9.7.23.final.bed

# Make bigbed file
/shen/shenlabstore3/ijones1/pre_PhD/copy_files_july_2019/home/ijones1/scripts/bedToBigBed  ${output}/Anchor.Merged.9.7.23.final.bed  $CHROM  ${output}/Anchor.Merged.9.7.23.final.bb
