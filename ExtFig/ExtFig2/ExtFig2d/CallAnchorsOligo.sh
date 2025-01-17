#/bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y
#$ -m ea
#$ -M Ian.Jones3@ucsf.edu

# Date: 4/13/22
# Goal: Combine short reads from NeuroHub project for determining Anchors regions

date
hostname

#### activate enviroment ####
source ~/.bashrc
conda activate deeptools

#### Var ####
# Location to output the files
output=/wynton/group/shen/NeuroHub.7.23.21/Final/PLAC/Anchor/Oligo


#### Downsample Oligo (3 reps) to a total of 60mil , 20mil per sample or 30 mil, 10 mil persample ####

#### Downsample IJ351 (56,558,630 reads, 28,279,315 paired reads) #### ####
echo 'starting on IJ351...'
# x='expr 20000000 / 56558630'
name=IJ351_YSXY14_100bp.fastp
x=42.3536 # ratio to get 10 million
y=42.7072 # ratio to get 20 million
bam=/wynton/group/shen/NeuroHub.7.23.21/PLAC-seq/feather_output/IJ351_YSXY14_100bp.fastp_20210927_102520/tempfiles

# Perform down sampling
samtools view -s $x -h -b ${bam}/${name}.shrt.bam > ${output}/${name}.10mil.shrt.bam
samtools sort ${output}/${name}.10mil.shrt.bam -o ${output}/${name}.10mil.shrt.sorted.bam

samtools view -s $y -h -b ${bam}/${name}.shrt.bam > ${output}/${name}.20mil.shrt.bam
samtools sort ${output}/${name}.20mil.shrt.bam -o ${output}/${name}.20mil.shrt.sorted.bam


#### Downsample IJ372 (81,366,086 reads, 40,683,043 paired reads) ####
echo 'starting on IJ372...'
# x='expr 20000000 / 81366086'
name=IJ372_YSJJ05_100bp.fastp
x=42.2458 # ratio to get 10 million
y=42.4916 # ratio to get 20 million
bam=/wynton/group/shen/NeuroHub.7.23.21/PLAC-seq/feather_output/IJ372_YSJJ05_100bp.fastp_20211103_224401/tempfiles

# Perform down sampling
samtools view -s $x -h -b ${bam}/${name}.shrt.bam > ${output}/${name}.10mil.shrt.bam
samtools sort ${output}/${name}.10mil.shrt.bam -o ${output}/${name}.10mil.shrt.sorted.bam

samtools view -s $y -h -b ${bam}/${name}.shrt.bam > ${output}/${name}.20mil.shrt.bam
samtools sort ${output}/${name}.20mil.shrt.bam -o ${output}/${name}.20mil.shrt.sorted.bam

#### Downsample IJ374 (55,739,362 reads, 27,869,681 paired reads) ####
echo 'starting on IJ372...'
# x='expr 20000000 / 81366086'
name=IJ374_YSJJ07_100bp.fastp
x=42.3588 # ratio to get 10 million
y=42.7176 # ratio to get 20 million
bam=/wynton/group/shen/NeuroHub.7.23.21/PLAC-seq/feather_output/IJ374_YSJJ07_100bp.fastp_20220112_192934/tempfiles

# Perform down sampling
samtools view -s $x -h -b ${bam}/${name}.shrt.bam > ${output}/${name}.10mil.shrt.bam
samtools sort ${output}/${name}.10mil.shrt.bam -o ${output}/${name}.10mil.shrt.sorted.bam

samtools view -s $y -h -b ${bam}/${name}.shrt.bam > ${output}/${name}.20mil.shrt.bam
samtools sort ${output}/${name}.20mil.shrt.bam -o ${output}/${name}.20mil.shrt.sorted.bam


####  Merge 10mil downsampled bam files ####
samtools merge ${output}/oligo.shrt.30mil.bam  ${output}/*.10mil.shrt.sorted.bam

samtools index ${output}/oligo.shrt.30mil.bam

bamCoverage -b ${output}/oligo.shrt.30mil.bam --normalizeUsing CPM -o ${output}/oligo.shrt.30mil.bw

#### Merge 20mil downsampled bam files ####
samtools merge ${output}/oligo.shrt.60mil.bam  ${output}/*.20mil.shrt.sorted.bam

samtools index ${output}/oligo.shrt.60mil.bam

bamCoverage -b ${output}/oligo.shrt.60mil.bam --normalizeUsing CPM -o ${output}/oligo.shrt.60mil.bw


### Call Peaks (based on Arima) ####
name30=oligo.shrt.30mil
name60=oligo.shrt.60mil
input_dir=${output}
genome=hs
macs_output=$output/peaks
CHROM=/shen/shenlabstore3/shared/reference_genome/hg38/hg38.chrom.sizes

# Call peaks
conda activate MACS
macs2 callpeak -t ${input_dir}/${name30}.bam -n ${name30} -g $genome --broad --nolambda --broad-cutoff 0.3 --outdir ${macs_output}/

macs2 callpeak -t ${input_dir}/${name60}.bam -n ${name60} -g $genome --broad --nolambda --broad-cutoff 0.3 --outdir ${macs_output}/

### Generate bigbed files ####
cut -f1-3 ${macs_output}/${name30}*.broadPeak > ${macs_output}/temp.bed
/shen/shenlabstore3/ijones1/pre_PhD/copy_files_july_2019/home/ijones1/scripts/bedToBigBed  ${macs_output}/temp.bed  $CHROM  ${macs_output}/${name30}.bb

cut -f1-3 ${macs_output}/${name60}*.broadPeak > ${macs_output}/temp.bed
/shen/shenlabstore3/ijones1/pre_PhD/copy_files_july_2019/home/ijones1/scripts/bedToBigBed  ${macs_output}/temp.bed  $CHROM  ${macs_output}/${name60}.bb
