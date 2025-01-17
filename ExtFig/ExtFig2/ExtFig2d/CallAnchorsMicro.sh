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
# Goal: Combine short reads from NeuroHub project for determining Anchors regions

date
hostname

#### activate enviroment ####
source ~/.bashrc
conda activate deeptools

#### Var ####
# Location to output the files
output=/wynton/group/shen/NeuroHub.7.23.21/Final/PLAC/Anchor/MG


#### Downsample vRG (3 reps) to a total of 30 mil, 10 mil persample ####

#### Downsample IJ361 #### ####
echo 'starting on IJ361...'
name=IJ361_YSJJ05_100bp.fastp
x=42.3143 # ratio to get 10 million
bam=/wynton/group/shen/NeuroHub.7.23.21/PLAC-seq/feather_output/IJ361_YSJJ05_100bp.fastp_20211103_222416/tempfiles

# Perform down sampling
# samtools view -s $x -h -b ${bam}/${name}.shrt.bam > ${output}/${name}.10mil.shrt.bam
# samtools sort ${output}/${name}.10mil.shrt.bam -o ${output}/${name}.10mil.shrt.sorted.bam


#### Downsample IJ375 ####
echo 'starting on IJ375...'
name=IJ375_YSJJ07_100bp.fastp
x=42.1888 # ratio to get 10 million
bam=/wynton/group/shen/NeuroHub.7.23.21/PLAC-seq/feather_output/IJ375_YSJJ07_100bp.fastp_20220112_192934/tempfiles

# Perform down sampling
# samtools view -s $x -h -b ${bam}/${name}.shrt.bam > ${output}/${name}.10mil.shrt.bam
# samtools sort ${output}/${name}.10mil.shrt.bam -o ${output}/${name}.10mil.shrt.sorted.bam

#### Downsample IJ394 ####
echo 'starting on IJ394...'
name=IJ394_YSJJ17_100bp.fastp
x=42.3818 # ratio to get 10 million
bam=/wynton/group/shen/NeuroHub.7.23.21/Final/PLAC/feather_output/IJ394_YSJJ17_100bp.fastp_20220714_113238/tempfiles

# Perform down sampling
samtools view -s $x -h -b ${bam}/${name}.shrt.bam > ${output}/${name}.10mil.shrt.bam
samtools sort ${output}/${name}.10mil.shrt.bam -o ${output}/${name}.10mil.shrt.sorted.bam




####  Merge 10mil downsampled bam files ####
samtools merge ${output}/MG.shrt.30mil.bam  ${output}/*.10mil.shrt.sorted.bam

samtools index ${output}/MG.shrt.30mil.bam

bamCoverage -b ${output}/MG.shrt.30mil.bam --normalizeUsing CPM -o ${output}/MG.shrt.30mil.bw



### Call Peaks (based on Arima) ####
name30=MG.shrt.30mil
input_dir=${output}
genome=hs
macs_output=$output/peaks
CHROM=/shen/shenlabstore3/shared/reference_genome/hg38/hg38.chrom.sizes

# Call peaks
conda activate MACS
macs2 callpeak -t ${input_dir}/${name30}.bam -n ${name30} -g $genome --broad --nolambda --broad-cutoff 0.3 --outdir ${macs_output}/

### Generate bigbed files ####
cut -f1-3 ${macs_output}/${name30}*.broadPeak > ${macs_output}/temp.bed
/shen/shenlabstore3/ijones1/pre_PhD/copy_files_july_2019/home/ijones1/scripts/bedToBigBed  ${macs_output}/temp.bed  $CHROM  ${macs_output}/${name30}.bb

