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

# Location of ATAC CPM bigwig
dir=/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/BIGWIG
# Location of peak bed files
dir2=/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison
# Location of methylation bigwig file
dir3=/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/Methylation/MG
# Cell type name
name=MG

computeMatrix reference-point -S ${dir}/Microglia_IJ310_JJ019_JJ076_hg38_CPM.bw ${dir3}/MG_CG.report.10.cov.bw -R ${dir2}/${name}.Both.peak.bed ${dir2}/${name}.ATAC.peak.bed ${dir2}/${name}.MR.peak.bed --referencePoint 'center' -b 1000 -a 1000 -bs 50 --sortRegions 'descend' --sortUsing 'mean' -o ${name}_Peaks.gz 

# plotHeatmap -m ${name}_Peaks.gz  -o ${name}_Peaks.pdf --colorMap Reds Blues --zMin 0 0  --zMax 2 1 --refPointLabel cCRE --regionsLabel 'LMAR' 'AR' 'LMR' --whatToShow 'plot' --colorList "#762A83,#9970AB,#C2A5CF"

plotProfile -m ${name}_Peaks.gz -out ${name}_Peaks.pdf --yMin 0 0  --yMax 2 1 --refPointLabel cCRE --regionsLabel 'LMAR' 'AR' 'LMR' --colors "#762A83" "#9970AB" "#C2A5CF"