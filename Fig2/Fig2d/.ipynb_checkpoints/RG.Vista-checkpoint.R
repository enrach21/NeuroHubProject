# Read in packages
library('bedtoolsr')
library('dplyr')
library('tidyr')
library(data.table)
library(ggplot2)
library("ComplexHeatmap")


# Read in Michael's Data
df <- read.table('/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/Scripts/ATAC-seq/Michael_tested/2022-10-21_for_Ian.tsv', header=T)
df1 <- df[,c('chr','start','end','element','cell','tissue')]
df1$outcome <- 'Pos'
df1$outcome[df1$tissue == 'NEGATIVE'] <- 'Neg'
df1


# Read in unique ATAC-seq peaks
ReadATAC <- function(dir, file){
    df<-read.table(paste0(dir,file))
    colnames(df) <- c('chr','start','end')
    return(df)
}

# Read in unique methylation peaks
ReadMRs <- function(dir, file){
    df<-read.table(paste0(dir,file), header=T)
    return(df)
}

# ATAC-seq Overlapping

# Location of called peaks
ATAC_dir <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/ATAC-seq/PEAKS/'

# Read in Data
MG_ATAC <- ReadATAC(ATAC_dir,'MG.overlap.optimal_peak.unique.bed')
nrow(MG_ATAC)
OPC_ATAC <- ReadATAC(ATAC_dir,'OPC.overlap.optimal_peak.unique.bed')
nrow(OPC_ATAC)
oRG_ATAC <- ReadATAC(ATAC_dir,'oRG.4reps.overlap.optimal_peak.unique.bed')
nrow(oRG_ATAC)
vRG_ATAC <- ReadATAC(ATAC_dir,'vRG.4reps.overlap.optimal_peak.unique.bed')
nrow(vRG_ATAC)



# Location of PLAC-seq files
DIR <- '/shen/shenlabstore3/ijones1/ShenAnalysis/NeuroHub/IntermediateFiles/PLAC_2023/bedpe/'

OC <- 'Oligo_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.bedpe'

MG <- 'Microglia_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.bedpe'

vRG <- 'vRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.bedpe'

oRG <- 'oRG_100bp.fastp.60mil.FRIP.2k.2.sig3Dinteractions.bedpe'

# Read in 
OC_DF <- read.table(paste0(DIR,OC), fill=T, header=T)

MG_DF <- read.table(paste0(DIR,MG), fill=T, header=T)

vRG_DF <- read.table(paste0(DIR,vRG), fill=T, header=T)

oRG_DF <- read.table(paste0(DIR,oRG), fill=T, header=T)


names <- c('MG','OC','oRG','vRG')
ATAC_list <- list(MG_ATAC,OPC_ATAC,oRG_ATAC,vRG_ATAC)
LMR_list <- list(MG_LMRs,OPC_LMRs,oRG_LMRs,vRG_LMRs)
PLAC_list <- list(MG_DF,OC_DF,oRG_DF,vRG_DF)

for (n in 1:length(names)) {
    temp <- bt.intersect(a=df1,b=ATAC_list[[n]], c=T)
    colnames(temp) <- c(colnames(df1),paste0(names[n],'_ATAC'))
           
    df1 <- temp
                        
    temp <- bt.intersect(a=df1,b=LMR_list[[n]], c=T)
    colnames(temp) <- c(colnames(df1),paste0(names[n],'_LMRs'))
                        
    df1 <- temp
    
    bed1 <- PLAC_list[[n]][1:3]
    bed2 <- PLAC_list[[n]][4:6]
    
    temp <- bt.intersect(a=df1,b=list(bed1,bed2), c=T)
    colnames(temp) <- c(colnames(df1),paste0(names[n],'_PLAC'))
    
    df1 <- temp
}

df2 <- df1[df1$cell=='RG',c(4,7:19)]
df2


# Re order
df2 <- df2[order(df2$outcome, df2$element),]
df2

# Get Data
mat <- as.matrix(df2[,3:14])
mat
# Change values to zero or 1
mat[,1:3][mat[,1:3] > 0] <- 'MG_Yes'
mat[,1:3][mat[,1:3] == 0] <- 'MG_No'

mat[,4:6][mat[,4:6] > 0] <- 'OC_Yes'
mat[,4:6][mat[,4:6] == 0] <- 'OC_No'

mat[,7:9][mat[,7:9] > 0] <- 'oRG_Yes'
mat[,7:9][mat[,7:9] == 0] <- 'oRG_No'

mat[,10:12][mat[,10:12] > 0] <- 'vRG_Yes'
mat[,10:12][mat[,10:12] == 0] <- 'vRG_No'

rownames(mat) <- df2$element
mat

colors = structure(c('#93AACF','#c9d4e7',
                    '#789678','#bbcabb',
                    '#E3803B','#f1bf9d',
                    '#6B4E30','#B5A697'), names = c("MG_Yes", "MG_No", "OC_Yes", "OC_No",
                                  "oRG_Yes",'oRG_No','vRG_Yes','vRG_No'))


anno <- HeatmapAnnotation(Vista = df2$outcome, 
            col = list(Vista=c("Pos"='red',"Neg"='blue')),
            gp = gpar(col = "white", lwd = 2),
            simple_anno_size = unit(.2, "cm"),
            show_annotation_name = c(Vista = FALSE))
anno

options(repr.plot.width=2.5, repr.plot.height=1.75)
pd <- Heatmap(t(mat), cluster_rows = FALSE, cluster_columns = FALSE,  
        row_split = rev(c('MG','MG','MG','OC','OC','OC','oRG','oRG','oRG','vRG','vRG','vRG')),
       row_labels = rep(c('ATAC','LMR','PLAC'), 4),
        show_heatmap_legend = FALSE, 
        col = colors, rect_gp = gpar(col = "white", lwd = .5),
       top_annotation = anno,
       row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7),
        row_title = NULL) 
pd

pdf('NH.Vista.Heatmap.pdf', height = 1.75, width = 2.5)
draw(pd)
dev.off()