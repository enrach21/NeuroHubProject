#### Create Control Files
library(gkmSVM)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
#check which genomes are installed
installed.genomes(splitNameParts=FALSE)
genome <- BSgenome.Hsapiens.UCSC.hg38.masked
genNullSeqs('/shen/shenlabstore3/ijones1/GKM_explain_test/vRG/Pos_90k/vRG.90000.peaks.bed', genome=genome, outputBedFN = '/shen/shenlabstore3/ijones1/GKM_explain_test/vRG/Neg_90k/vRG.neg.peaks.bed', outputPosFastaFN = '/shen/shenlabstore3/ijones1/GKM_explain_test/vRG/Pos_90k/vRG.90000.peaks.fa', outputNegFastaFN = '/shen/shenlabstore3/ijones1/GKM_explain_test/vRG/Neg_90k/vRG.neg.peaks.fa', nMaxTrials=20, xfold=2)