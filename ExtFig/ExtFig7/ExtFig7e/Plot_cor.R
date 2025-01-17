library('bedtoolsr')
library('dplyr')
library(ggplot2)
library("UpSetR")
library("reshape2")
library('ggpointdensity')
library(viridis)

# Read in model
df <- read.csv('MG.AD.8.14.24.csv', row.names = 1)

 cor.test(df$explain_score, df$delta_score, method = "pearson", conf.level = 0.95)
ggplot(df, aes(x=explain_score,y=delta_score)) +
    # geom_point()
      geom_pointdensity(size = .2) +
      scale_color_viridis() +
      theme_bw() +
       theme( plot.title = element_text(color="Black", size=8, hjust = 0.5),
       # axis.text.x = element_text( angle=45, hjust = 1),
       axis.text = element_text(size = 8),
       legend.position = "none",
       axis.title=element_text(size=8)) +
    xlab('GKM Explain') +
    ylab("GKM Delta") + geom_smooth(method=lm, se=FALSE)
ggsave('MG.explain.delta.pdf', height = 1.25, width =1.25)

 cor.test(df$ISM_score, df$delta_score, method = "pearson", conf.level = 0.95)
ggplot(df, aes(x=ISM_score,y=delta_score)) +
    # geom_point()
      geom_pointdensity(size = .2) +
      scale_color_viridis() +
      theme_bw() +
       theme( plot.title = element_text(color="Black", size=8, hjust = 0.5),
       # axis.text.x = element_text( angle=45, hjust = 1),
       axis.text = element_text(size = 8),
       legend.position = "none",
       axis.title=element_text(size=8)) +
    xlab('ISM') +
    ylab("Delta") + geom_smooth(method=lm, se=FALSE)
ggsave('MG.ISM.delta.pdf', height = 1.25, width =1.25)

 cor.test(df$explain_score, df$ISM_score, method = "pearson", conf.level = 0.95)
ggplot(df, aes(x=ISM_score,y=explain_score)) +
    # geom_point()
      geom_pointdensity(size = .2) +
      scale_color_viridis() +
      theme_bw() +
       theme( plot.title = element_text(color="Black", size=8, hjust = 0.5),
       # axis.text.x = element_text( angle=45, hjust = 1),
       axis.text = element_text(size = 8),
       legend.position = "none",
       axis.title=element_text(size=8)) +
    xlab('ISM') +
    ylab("Explain") + geom_smooth(method=lm, se=FALSE)
ggsave('MG.explain.ISM.pdf', height = 1.25, width =1.25)