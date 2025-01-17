library(PRROC)

#### vRG ####

pos <- read.table('vRG_Pos_chr2_predict_output.txt')
dim(pos)
head(pos)

neg <- read.table('vRG_Neg_chr2_predict_output.txt')
dim(neg)
head(neg)

pr <- pr.curve(pos$V2, neg$V2, curve=T)
plot(pr)
roc <- roc.curve(pos$V2, neg$V2, curve=T)
plot(roc)

pdf('vRG_AUCROC.pdf', height = 4, width = 3.5)
plot(roc, legend=FALSE, color='#2166AC', cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5)
dev.off()

plot(roc, legend=FALSE, color='#2166AC',     cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)