#### Get target Genes of motifs ####

### The following done on the server at
cd /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/CytoScapePlot/LHX2

# oRG
# PLAC
bedtools intersect -a oRG_LHX2_motif.bed -b /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.V2.bedpe -wa -wb | cut -f 1-5,22 > oRG_LHX2_motif_PLAC.bed

# TSS
bedtools intersect -a oRG_LHX2_motif.bed -b /shen/shenlabstore3/shared/PLAC-seq_analysis/utils/TssFiles/gencode.v38.tss.1000bp.update.bed -wa -wb | cut -f 1-5,12 > oRG_LHX2_motif_TSS.bed


# vRG
# PLAC
bedtools intersect -a vRG_LHX2_motif.bed -b /shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.mod.V2.bedpe -wa -wb | cut -f 1-5,22 > vRG_LHX2_motif_PLAC.bed

# TSS
bedtools intersect -a vRG_LHX2_motif.bed -b /shen/shenlabstore3/shared/PLAC-seq_analysis/utils/TssFiles/gencode.v38.tss.1000bp.update.bed -wa -wb | cut -f 1-5,12 > vRG_LHX2_motif_TSS.bed
