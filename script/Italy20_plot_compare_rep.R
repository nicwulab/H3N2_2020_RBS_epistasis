#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
require(cowplot)

PlotCompareFit_Rep <- function(Italy20_data, graphname){
    textsize <- 7
    colorscale  <- c(brewer.pal(9,"GnBu"))
    p <- ggplot() +
           geom_rect(data=NULL,aes(xmin=0,xmax=Inf,ymin=0,ymax=Inf), color=NA, fill='grey85', alpha=1) +
           geom_point(data=Italy20_data, aes(x=log10(Rep1Enrich), y=log10(Rep2Enrich),fill=mutclass),
                      pch=21, size=1, alpha=1, stroke=0.2) +
           geom_vline(xintercept=0, lty="11", col = 'black') +
           geom_hline(yintercept=0, lty="11", col = 'black') +
           scale_fill_gradient2(limits=c(1,10),midpoint=4.5,low=colorscale[1],mid=colorscale[5],high=colorscale[9],
                            breaks=c(1,10),labels=c(1,10),
                            name='# of substitutions',
                            guide=guide_colorbar(title.position="top",
                                                 title.vjust=0,
                                                 label.position="bottom")) +
           theme_cowplot(12) +
           theme(axis.title=element_text(size=textsize,face="bold",family="Arial"),
                 axis.text=element_text(size=textsize,face="bold",family="Arial"),
                 legend.title=element_blank(),
                 legend.key.size=unit(0.1,'in'),
                 legend.spacing.x=unit(0.03, 'in'),
                 legend.text=element_text(size=textsize,face="bold",family="Arial"),
                 legend.position='top') +
           labs(x=expression(bold(log['10']~enrich~'(Rep 1)')),y=expression(bold(log['10']~enrich~'(Rep 2)')))
    ggsave(graphname,p,height=2.8,width=2.5,dpi=600, bg='white')
    } 

Italy20_data <- read_tsv('results/Ita20HA_MultiMutLib_filtered.tsv') %>%
                  filter(grepl('L194P', mut)) %>%
                  mutate(mutclass=mutclass-1)
print (Italy20_data)
PlotCompareFit_Rep(Italy20_data,'graph/Italy20_mutlib_rep_compare.png')
