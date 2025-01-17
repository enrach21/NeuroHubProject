#!/bin/bash
#$ -l h_rt=72:0:0
#$ -l mem_free=16G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -M Ian.Jones3@ucsf.edu

source ~/.bashrc
conda activate MAPS_update

#######################################
###INPUT VARIABLES
#######################################
Rscript_path="/shen/shenlabstore3/ijones1/dependencies/conda/envs/encode.test/envs/MAPS_update/bin/Rscript"
dir_name="/wynton/group/shen/NeuroHub.7.23.21/PLAC-seq/HPREP/2023_10mil_usable/HPRep_output/stage3/"
tune_sample_1="IJ361_YSJJ05_100bp.fastp.norm.5k.normalized.txt"
tune_sample_2="IJ386_YSJJ14_100bp.fastp.norm.5k.normalized.txt"
chr_count="22"
bin_size="5000"
binning_range="2000000"
output_name="NeuroHub.10million.9_8_23.Norm"
seed="1"
#######################################
# Changed it to point at HPRep
$Rscript_path /shen/shenlabstore3/ijones1/dependencies/HPRep/bin/HPRep/HPRep.R $dir_name $tune_sample_1 $tune_sample_2 $chr_count $bin_size $binning_range $output_name $seed




