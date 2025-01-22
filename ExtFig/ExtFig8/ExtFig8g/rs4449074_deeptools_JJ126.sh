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
bamCoverage -b JJ126_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bam --binSize 1 --scaleFactor 0.9263 -r chr2:199258295:199264414 -o rs4449074.JJ126_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bw
fi

if [ ${SGE_TASK_ID} == 2 ]; then 
echo 'Starting 2'
bamCoverage -b JJ126_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bam --binSize 1 --normalizeUsing CPM -r chr2:199258295:199264414 -o rs4449074.CPM.JJ126_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bw
fi

if [ ${SGE_TASK_ID} == 3 ]; then 
echo 'Starting 3' 
bamCoverage -b JJ126_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bam --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2150570000 --ignoreForNormalization chrX -r chr2:199258295:199264414 -o rs4449074.RPGC.JJ126_YSJJ17.fastp.R1.trim.merged.srt.nodup.no_chrM_MT.bw
fi



