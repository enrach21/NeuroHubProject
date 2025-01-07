#/bin/bash
#$ -l h_rt=1:0:0
#$ -l mem_free=4G
#$  -cwd
#$ -l scratch=100G
#$ -j y
#$ -r y

date
hostname

# activate enviroment:
source ~/.bashrc

awk -F'\t' -v OFS='\t' 'NR { $(NF+1)=NR} 1' MG.overlap.optimal_peak.unique.bed > MG.overlap.optimal_peak.name.unique.bed

awk -F'\t' -v OFS='\t' 'NR { $(NF+1)=NR} 1' OPC.overlap.optimal_peak.unique.bed > OPC.overlap.optimal_peak.name.unique.bed


/shen/shenlabstore3/ijones1/dependencies/MAPS_update/MAPS-master/bin/utils/genomic_features_generator/bigWigAverageOverBed /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/BIGWIG/Microglia_IJ310_JJ019_JJ076_hg38_CPM.bw /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/PEAKS/Merged.4.11.23.overlap.optimal_peak.final2.bed MG.11.21.23.overlap.optimal_peak.final.txt -bedOut=MG.11.21.23.overlap.optimal_peak.count.bed 

/shen/shenlabstore3/ijones1/dependencies/MAPS_update/MAPS-master/bin/utils/genomic_features_generator/bigWigAverageOverBed /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/BIGWIG/Oligo_IJ309_JJ020_JJ077_hg38_CPM.bw /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/PEAKS/Merged.4.11.23.overlap.optimal_peak.final2.bed OPC.11.21.23.overlap.optimal_peak.final.txt -bedOut=OPC.11.21.23.overlap.optimal_peak.count.bed 


/shen/shenlabstore3/ijones1/dependencies/MAPS_update/MAPS-master/bin/utils/genomic_features_generator/bigWigAverageOverBed /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/BIGWIG/vRG_IJ308_IJ389_IJ393_JJ126_hg38_CPM.bw /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/PEAKS/Merged.4.11.23.overlap.optimal_peak.final2.bed vRG.4reps.8.15.23.overlap.optimal_peak.final.txt -bedOut=vRG.4reps.8.15.23.overlap.optimal_peak.count.bed 

/shen/shenlabstore3/ijones1/dependencies/MAPS_update/MAPS-master/bin/utils/genomic_features_generator/bigWigAverageOverBed /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/BIGWIG/oRG_IJ307_IJ388_IJ392_JJ125_hg38_CPM.bw /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/PEAKS/Merged.4.11.23.overlap.optimal_peak.final2.bed oRG.4reps.8.15.23.overlap.optimal_peak.final.txt -bedOut=oRG.4reps.8.15.23.overlap.optimal_peak.count.bed 