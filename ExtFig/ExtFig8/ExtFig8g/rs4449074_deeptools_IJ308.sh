#/bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=16G
#$  -cwd
#$ -l scratch=500G
#$ -j y
#$ -r y
#$ -t 1-3

date
hostname

# activate enviroment:
source ~/.bashrc
conda activate deeptools


if [ ${SGE_TASK_ID} == 1 ]; then 
echo 'Starting 1'
bamCoverage -b IJ308_YSXY12.fastp.R1.trim.srt.nodup.no_chrM_MT.bam --binSize 1 --scaleFactor 2.1061 -r chr2:199258295:199264414 -o rs4449074.IJ308_YSXY12.fastp.R1.trim.srt.nodup.no_chrM_MT.bw
fi

if [ ${SGE_TASK_ID} == 2 ]; then 
echo 'Starting 2'
bamCoverage -b IJ308_YSXY12.fastp.R1.trim.srt.nodup.no_chrM_MT.bam --binSize 1 --normalizeUsing CPM -r chr2:199258295:199264414 -o rs4449074.CPM.IJ308_YSXY12.fastp.R1.trim.srt.nodup.no_chrM_MT.bw
fi

if [ ${SGE_TASK_ID} == 3 ]; then 
echo 'Starting 3' 
bamCoverage -b IJ308_YSXY12.fastp.R1.trim.srt.nodup.no_chrM_MT.bam --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2150570000 --ignoreForNormalization chrX -r chr2:199258295:199264414 -o rs4449074.RPGC.IJ308_YSXY12.fastp.R1.trim.srt.nodup.no_chrM_MT.bw
fi



