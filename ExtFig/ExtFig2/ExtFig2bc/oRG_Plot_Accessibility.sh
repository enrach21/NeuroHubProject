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

dir=/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/BIGWIG
dir2=/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison
# Location of methylation bigwig file
dir3=/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/Methylation/oRG
name=oRG

computeMatrix reference-point -S ${dir}/oRG_IJ307_IJ388_IJ392_JJ125_hg38_CPM.bw ${dir3}/oRG_CG.report.10.cov.bw -R ${dir2}/${name}.Both.peak.bed ${dir2}/${name}.ATAC.peak.bed ${dir2}/${name}.MR.peak.bed --referencePoint 'center' -b 1000 -a 1000 -bs 50 --sortRegions 'descend' --sortUsing 'mean' -o ${name}_Peaks.gz 

# plotHeatmap -m ${name}_Peaks.gz  -o ${name}_Peaks.pdf --colorMap RdBu --refPointLabel cCRE --regionsLabel 'LMAR' 'AR' 'LMR'


plotProfile -m ${name}_Peaks.gz -out ${name}_Peaks.pdf --yMin 0 0  --yMax 2 1 --refPointLabel cCRE --regionsLabel 'LMAR' 'AR' 'LMR' --colors "#B2182B" "#D6604D" "#F4A582"