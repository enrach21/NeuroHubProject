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
output=/wynton/group/shen/NeuroHub.7.23.21/Final/PLAC/Anchor/oRG


#### Downsample oRG (3 reps) to a total of 30 mil, 10 mil persample ####

#### Downsample IJ356  #### ####
echo 'starting on IJ356...'
name=IJ356_YSXY14_100bp.fastp
x=42.9832228748 # ratio to get 10 million
bam=/wynton/group/shen/NeuroHub.7.23.21/PLAC-seq/feather_output/IJ356_YSXY14_100bp.fastp_20210927_103218/tempfiles

# Perform down sampling
samtools view -s $x -h -b ${bam}/${name}.shrt.bam > ${output}/${name}.10mil.shrt.bam
samtools sort ${output}/${name}.10mil.shrt.bam -o ${output}/${name}.10mil.shrt.sorted.bam


#### Downsample IJ363 ####
echo 'starting on IJ363...'
name=IJ363_YSJJ05_100bp.fastp
x=42.3409 # ratio to get 10 million
bam=/wynton/group/shen/NeuroHub.7.23.21/PLAC-seq/feather_output/IJ363_YSJJ05_100bp.fastp_20211103_222417/tempfiles

# Perform down sampling
samtools view -s $x -h -b ${bam}/${name}.shrt.bam > ${output}/${name}.10mil.shrt.bam
samtools sort ${output}/${name}.10mil.shrt.bam -o ${output}/${name}.10mil.shrt.sorted.bam


#### Downsample IJ385 ####
echo 'starting on IJ385...'
name=IJ385_YSJJ14_100bp.fastp
x=42.4947 # ratio to get 10 million
bam=/wynton/group/shen/NeuroHub.7.23.21/PLAC-seq/feather_output/IJ385_YSJJ14_100bp.fastp_20220414_154619/tempfiles

# Perform down sampling
samtools view -s $x -h -b ${bam}/${name}.shrt.bam > ${output}/${name}.10mil.shrt.bam
samtools sort ${output}/${name}.10mil.shrt.bam -o ${output}/${name}.10mil.shrt.sorted.bam



####  Merge 10mil downsampled bam files ####
samtools merge ${output}/oRG.shrt.30mil.bam  ${output}/*.10mil.shrt.sorted.bam

samtools index ${output}/oRG.shrt.30mil.bam

bamCoverage -b ${output}/oRG.shrt.30mil.bam --normalizeUsing CPM -o ${output}/oRG.shrt.30mil.bw



### Call Peaks (based on Arima) ####
name30=oRG.shrt.30mil
name60=oRG.shrt.60mil
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

