#### Packages

library('bedtoolsr')
library('dplyr')
library('tidyr')
library(ggplot2)
library("UpSetR")
library("reshape2")


# Up in Ref
vRG <- read.table('/shen/shenlabstore3/ijones1/GKM_explain_test/vRG/AD_scores_update/vRG.AD.Sig.Ref.txt')$x

oRG <- read.table('/shen/shenlabstore3/ijones1/GKM_explain_test/oRG/AD_scores_update/oRG.AD.Sig.Ref.txt')$x

OPC <- read.table('/shen/shenlabstore3/ijones1/GKM_explain_test/OPC/AD_scores_update/OPC.AD.Sig.Ref.txt')$x

MG <- read.table('/shen/shenlabstore3/ijones1/GKM_explain_test/MG/AD_scores_update/MG.AD.Sig.Ref.txt')$x

listInput <- list(vRG = vRG , oRG = oRG , OPC=OPC, MG=MG)

pdf('AD.Upset.Ref.8.27.24.pdf', height=3, width=3.5)
upset(fromList(listInput), order.by = "freq")
dev.off()


# up in alt

vRG <- read.table('/shen/shenlabstore3/ijones1/GKM_explain_test/vRG/AD_scores_update/vRG.AD.Sig.Alt.txt')$x

oRG <- read.table('/shen/shenlabstore3/ijones1/GKM_explain_test/oRG/AD_scores_update/oRG.AD.Sig.Alt.txt')$x

OPC <- read.table('/shen/shenlabstore3/ijones1/GKM_explain_test/OPC/AD_scores_update/OPC.AD.Sig.Alt.txt')$x

MG <- read.table('/shen/shenlabstore3/ijones1/GKM_explain_test/MG/AD_scores_update/MG.AD.Sig.Alt.txt')$x

listInput <- list(vRG = vRG , oRG = oRG , OPC=OPC, MG=MG)


pdf('AD.Upset.Alt.8.27.24.pdf', height=3, width=3.5)
upset(fromList(listInput), order.by = "freq")
dev.off()