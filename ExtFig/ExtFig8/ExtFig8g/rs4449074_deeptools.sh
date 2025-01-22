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



# samtools view JJ126_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bam "chr2:199261161-199261604" -b > rs4449074.JJ126_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bam

bamCoverage -b rs4449074.JJ126_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bam --binSize 1 --scaleFactor 0.9263 -o rs4449074.JJ126_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bw

# samtools view IJ393_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bam "chr2:199261161-199261604" -b > rs4449074.IJ393_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bam

bamCoverage -b rs4449074.IJ393_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bam --binSize 1 --scaleFactor 0.7402 -o rs4449074.IJ393_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bw

# samtools view IJ389_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bam "chr2:199261161-199261604" -b > rs4449074.IJ389_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bam

bamCoverage -b rs4449074.IJ389_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bam --binSize 1 --scaleFactor 0.6800 -o rs4449074.IJ389_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bw

# samtools view IJ308_YSXY12.fastp.R1.trim.srt.nodup.no_chrM_MT.bam "chr2:199261161-199261604" -b > rs4449074.IJ308_YSXY12.fastp.R1.trim.srt.nodup.no_chrM_MT.bam

bamCoverage -b rs4449074.IJ308_YSXY12.fastp.R1.trim.srt.nodup.no_chrM_MT.bam --binSize 1 --scaleFactor 2.1061 -o rs4449074.IJ308_YSXY12.fastp.R1.trim.srt.nodup.no_chrM_MT.bw



