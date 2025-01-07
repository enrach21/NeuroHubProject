library('bedtoolsr')
library('dplyr')
library('tidyr')
library(data.table)
library(ggplot2)
library(gprofiler2)
library("ggrepel")
library(forcats)

MG <- read.table('MG.genes.txt')
OPC <- read.table('OPC.genes.txt')
vRG <- read.table('vRG.genes.txt')
oRG <- read.table('oRG.genes.txt')
genes <- unique(c(MG$V1, OPC$V1, vRG$V1, oRG$V1))
genes <- genes[!is.na(genes)]
length(genes)

g.term3 <- gost(query = list("VISTA"=genes), organism = "hsapiens", ordered_query = F, significant = T, correction_method = "g_SCS", sources = "GO:BP", evcodes = TRUE)

res <- g.term3$result

fwrite(res, 'Vista.Targets.go.csv')