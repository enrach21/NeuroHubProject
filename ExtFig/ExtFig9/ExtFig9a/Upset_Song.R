
source ~/.bashrc
conda activate intervene


# Use Intervene to obtain groups of cCREs
intervene upset -i IPC*.uniq.bed eN*.uniq.bed iN*.uniq.bed /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/vRG.Both.peak.bed /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/oRG.Both.peak.bed \
                    --output ./Intervene_Song_NH_V3 --save-overlaps 
                    
intervene venn -i IPC*.uniq.bed eN*.uniq.bed iN*.uniq.bed /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/vRG.Both.peak.bed /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/oRG.Both.peak.bed \
                    --output ./Intervene_Venn_Song_NH_V3 --save-overlaps 


#### Create bed file to sample from
cat IPC*.uniq.bed eN*.uniq.bed iN*.uniq.bed /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/vRG.Both.peak.bed /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq_MR_comparison/oRG.Both.peak.bed | bedtools sort | bedtools merge > Song.2020.NH.V3.ATAC.bed