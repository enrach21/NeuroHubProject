# activate enviroment:
source ~/.bashrc
conda activate homer

# Names
IN_DIR=/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq/DARs/vRG_oRG
BED=vRG_peaks.vRG_oRG.diffbind_analysis.bed
BED2=oRG_peaks.vRG_oRG.diffbind_analysis.bed
OUTPUT=/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq/DARs/vRG_oRG/HOMER
DIR=vRG_genome.Background_10_5_23
DIR2=oRG_genome.Background_10_5_23

# Commands

mkdir ${OUTPUT}/${DIR}
mkdir ${OUTPUT}/${DIR2}


findMotifsGenome.pl ${IN_DIR}/${BED} hg38 ${OUTPUT}/${DIR} -size given 


findMotifsGenome.pl ${IN_DIR}/${BED2} hg38 ${OUTPUT}/${DIR2} -size given 