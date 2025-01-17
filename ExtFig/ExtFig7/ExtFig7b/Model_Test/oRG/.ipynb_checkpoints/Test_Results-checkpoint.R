library(PRROC)

#### oRG ####
pos <- read.table('oRG_Pos_chr2_predict_output.txt')
dim(pos)
head(pos)

neg <- read.table('oRG_Neg_chr2_predict_output.txt')
dim(neg)
head(neg)

pr <- pr.curve(pos$V2, neg$V2, curve=T)
plot(pr)
roc <- roc.curve(pos$V2, neg$V2, curve=T)
plot(roc)

pdf('oRG_AUCROC.pdf', height = 4, width = 3.5)
plot(roc, legend=FALSE, color='#B2182B', cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5)
dev.off()