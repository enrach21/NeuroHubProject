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


multiBamSummary BED-file --BED /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/PEAKS/Merged.4.11.23.overlap.optimal_peak.final2.bed --bamfiles ${dir}/*.no_chrM_MT.bam -o NH.Merged.9.12.23.npz --scalingFactors NH.Merged.9.12.23.txt --outRawCounts NH.Merged.9.12.23.tab

plotCorrelation \
    -in NH.Merged.9.12.23.npz \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap coolwarm \
    -o NH.Merged.9.12.23.heatmap_SpearmanCorr_readCounts.pdf --plotNumbers \
    --outFileCorMatrix NH.Merged.9.12.2.3SpearmanCorr_readCounts.tab \
    --labels oRG_R1 vRG_R1 OPC_R1 MG_R1 oRG_V2 vRG_R2 oRG_R3 vRG_R3 MG_R2 OPC_R2 MG_R3 OC_R3 oRG_R4 vRG_R4 \
    --plotHeight 5 \
    --plotWidth 5