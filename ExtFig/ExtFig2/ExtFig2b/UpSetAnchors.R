# Date: 9/8/23
# Author: Ian Jones
# Goal: Make an upset plot of Anchors


#### Packages ####
library('bedtoolsr')
library('dplyr')
library(data.table)
library(ggplot2)
library("UpSetR")
library("ComplexHeatmap")
library(gprofiler2)

#### Functions ####

get_peak <- function(file, bin){
    temp <- read.table(file)[,1:3]
    binDF <- read.table(bin)
    temp2<-bt.intersect(a=binDF, b=temp,wa=T,u=T)
    peaks <- paste0(temp2$V1,':',temp2$V2,'-',temp2$V3)
    return(peaks)
}


#### Read in meta data
df<- read.csv('AnchorTable2KB_0.01.tsv')


if (nrow(df)==4) {
    df$peak <- NA
    listInput <- c()
    for (x in 1:nrow(df)) { 
        listInput <- c(listInput, setNames(list(get_peak(df$File[x], df$Bin[x])), df$Name[x]))
    }
}


pdf('NH.Anchors.9.8.23.pdf', width=8, height = 6)
upset(fromList(listInput), main.bar.color = "grey45",matrix.color="grey45"
      , order.by = "freq", text.scale=1.3,
      sets.x.label = 'Anchor Bins', mainbar.y.label='Anchor Bin Intersections',
     queries = list(
         list(query = intersects, 
    params = list(df$Name[1]), color = df$Color[1], active = T), 
        list(query = intersects, 
    params = list(df$Name[2]), color = df$Color[2], active = T),
     list(query = intersects, 
    params = list(df$Name[3]), color = df$Color[3], active = T),
     list(query = intersects, 
    params = list(df$Name[4]), color = df$Color[4], active = T)))
dev.off()

