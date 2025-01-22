#### Packages ####
library('bedtoolsr')
library('dplyr')
library('tidyr')
library(ggplot2)
library("UpSetR")
library("reshape2")


#### For each cell type ####
df <- read.csv('vRG.SCZ.8.12.24.csv')[,3:65]

write.table(unique(df[df$ISM_pval < 0.05 & df$explain_pval < 0.05 & df$delta_p_val < 0.05 
                & df$observed_active_allele == 'effect' 
                & df$vRG_ATAC > 0,'name']),'SCZ.vRG.Ref.txt')

write.table(unique(df[df$ISM_pval < 0.05 & df$explain_pval < 0.05 & df$delta_p_val < 0.05 
                & df$observed_active_allele == 'noneffect' 
                & df$vRG_ATAC > 0,'name']),'SCZ.vRG.Alt.txt')

# Plot Ref higher

vRG <- read.table('/shen/shenlabstore3/ijones1/GKM_explain_test/vRG/SCZ_score_update/SCZ.vRG.Ref.txt')$x
length(vRG)
oRG <- read.table('/shen/shenlabstore3/ijones1/GKM_explain_test/oRG/SCZ_scores_update/SCZ.oRG.Ref.txt')$x
length(oRG)
OPC <- read.table('/shen/shenlabstore3/ijones1/GKM_explain_test/OPC/SCZ_scores_update/SCZ.OPC.Ref.txt')$x
length(OPC)
MG <- read.table('/shen/shenlabstore3/ijones1/GKM_explain_test/MG/SCZ_scores_update/SCZ.MG.Ref.txt')$x
length(MG)

listInput <- list(vRG = vRG , oRG = oRG , OPC=OPC, MG=MG)

pdf('SCZ.Upset.Ref.8.27.24.pdf', height=3, width=3.5)
upset(fromList(listInput), order.by = "freq")
dev.off()

# Plot alt higher
vRG <- read.table('/shen/shenlabstore3/ijones1/GKM_explain_test/vRG/SCZ_score_update/SCZ.vRG.Alt.txt')$x

oRG <- read.table('/shen/shenlabstore3/ijones1/GKM_explain_test/oRG/SCZ_scores_update/SCZ.oRG.Alt.txt')$x

OPC <- read.table('/shen/shenlabstore3/ijones1/GKM_explain_test/OPC/SCZ_scores_update/SCZ.OPC.Alt.txt')$x

MG <- read.table('/shen/shenlabstore3/ijones1/GKM_explain_test/MG/SCZ_scores_update/SCZ.MG.Alt.txt')$x

listInput <- list(vRG = vRG , oRG = oRG , OPC=OPC, MG=MG)

pdf('SCZ.Upset.Alt.8.27.24.pdf', height=3, width=3.5)
upset(fromList(listInput), order.by = "freq")
dev.off()
