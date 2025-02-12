library('bedtoolsr')
library('dplyr')
library('tidyr')
library('ggplot2')
library(data.table)
library(readxl)

df <- read.csv('Luciferase.Pilot.11.6.24 - Chimp_Var.csv', row.names = 1, header=T)
df2 <- df[c('T1','T2','T3'),-2]
df2

df3 <- melt(df2)

  # grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

df4 <- data_summary(df3, varname="value", 
                    groupnames='variable')
# Convert dose to a factor variable
head(df4)

options(repr.plot.width = 4, repr.plot.height = 3)
p <- ggplot(df3, aes(x=variable, y=value, fill=variable)) +
    stat_summary( fun.data = mean_se, geom = "col", color ="black", position = position_dodge(1)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(1)) +
      geom_jitter(aes(),
              size = 1, shape=1) +
  theme_classic() +
  scale_fill_manual(values = c("grey", rep(c("#B2182B", '#FDDBC7'),3))) +
 guides(fill="none") +
 xlab('') + ylab('')
p
ggsave('Chimp.pdf',width=4, height = 3)