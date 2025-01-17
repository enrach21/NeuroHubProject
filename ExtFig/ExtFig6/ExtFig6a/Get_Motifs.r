#### Get motifs ####
# Read in Motifs
df <- read_excel('LHX2_HUMAN.H11MO.0.A_LHX2_HUMAN.H11MO.0.A_overview.xlsx') 

# Filter for motifs bound in oRG and >1 log2FC
oRG_motif <- df[df$oRG_bound==1 & df$oRG_vRG_log2fc >=1 ,]
oRG_motif <- oRG_motif[,c('TFBS_chr','TFBS_start','TFBS_end','oRG_score','oRG_vRG_log2fc')]
dim(oRG_motif)

write.table(oRG_motif, 'oRG_LHX2_motif.bed',col.names = F, row.names = F, quote = F, sep= '\t')

vRG_motif <- df[df$vRG_bound==1 & df$oRG_vRG_log2fc <= -1 ,]
vRG_motif <- vRG_motif[,c('TFBS_chr','TFBS_start','TFBS_end','vRG_score','oRG_vRG_log2fc')]

write.table(vRG_motif, 'vRG_LHX2_motif.bed',col.names = F, row.names = F, quote = F, sep= '\t')

# Read in Motifs
df <- read_excel('ASCL1_HUMAN.H11MO.0.A_ASCL1_HUMAN.H11MO.0.A_overview.xlsx') 

# Filter for motifs bound in oRG and >1 log2FC
oRG_motif <- df[df$oRG_bound==1 & df$oRG_vRG_log2fc >=1 ,]
oRG_motif <- oRG_motif[,c('TFBS_chr','TFBS_start','TFBS_end','oRG_score','oRG_vRG_log2fc')]

write.table(oRG_motif, 'oRG_ASCL1_motif.bed',col.names = F, row.names = F, quote = F, sep= '\t')

vRG_motif <- df[df$vRG_bound==1 & df$oRG_vRG_log2fc <= -1 ,]
vRG_motif <- vRG_motif[,c('TFBS_chr','TFBS_start','TFBS_end','vRG_score','oRG_vRG_log2fc')]
dim(vRG_motif)
write.table(vRG_motif, 'vRG_ASCL1_motif.bed',col.names = F, row.names = F, quote = F, sep= '\t')
