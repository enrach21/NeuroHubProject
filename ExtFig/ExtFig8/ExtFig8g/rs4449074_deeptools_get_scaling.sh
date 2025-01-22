#/bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y

date
hostname

# activate enviroment:
source ~/.bashrc
conda activate deeptools

dir=/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/BAM

multiBamSummary BED-file --BED /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/PEAKS/vRG.4reps.overlap.optimal_peak.unique.bed --bamfiles ${dir}/IJ308_YSXY12.fastp.R1.trim.srt.nodup.no_chrM_MT.bam ${dir}/IJ389_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bam ${dir}/IJ393_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bam ${dir}/JJ126_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bam -o NH.vRG.6.28.24.npz --scalingFactors NH.scale.vRG.6.28.24.txt --outRawCounts NH.vRG.6.28.24.tab

