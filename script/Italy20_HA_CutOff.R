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

filtered_data <- read.table(file = '../results/Ita20HA_MultiMutLib_filtered.tsv', sep = '\t', header = TRUE)
rep1_filtered_freq <- vector()
counter <-0
for (i in filtered_data$Rep1Count)
  {counter <- counter +1
  rep1_filtered_freq[counter]<-i/sum(filtered_data$Rep1Count)}
filtered_data$Rep1FilteredFreq <- rep1_filtered_freq
counter <-0
rep2_filtered_freq <- vector()
for (i in filtered_data$Rep2Count)
  {counter <- counter +1
  rep2_filtered_freq[counter]<-i/sum(filtered_data$Rep2Count)}
filtered_data$Rep2FilteredFreq <- rep2_filtered_freq

enrichment_correlation <-function (df, freq_cutoff_exponent){
  dataset_cutoff <- df[which(df$Rep1FilteredFreq > 10^freq_cutoff_exponent & df$Rep2FilteredFreq > 10^freq_cutoff_exponent),]
  dataset_after_cutoff<- dataset_cutoff%>% filter(grepl('L194P', mut))
  post_cutoff_input_sum <- sum(dataset_after_cutoff$InputCount)
  post_cutoff_rep1_sum <- sum(dataset_after_cutoff$Rep1Count)
  post_cutoff_rep2_sum <- sum(dataset_after_cutoff$Rep2Count)
  input_freq <- dataset_after_cutoff$InputCount/post_cutoff_input_sum
  rep1_freq <- dataset_after_cutoff$Rep1Count/post_cutoff_rep1_sum
  rep2_freq <- dataset_after_cutoff$Rep2Count/post_cutoff_rep2_sum
  rep1_enrichment <- rep1_freq/input_freq
  rep2_enrichment <- rep2_freq/input_freq
  enrichment_df <- data.frame(dataset_after_cutoff$mut, dataset_after_cutoff$mutclass, rep1_freq, rep2_freq, rep1_enrichment, rep2_enrichment)
  colnames(enrichment_df)[1] <- "mutants_containing_L194P"
  colnames(enrichment_df)[2] <- "mutation_counts"
  return(enrichment_df)
  }

graph_correlation <- function (enrichment_df){
  print(paste("Pearson Cor: ", round(cor(enrichment_df$rep1_enrichment, enrichment_df$rep2_enrichment),4), sep=''))
  textsize <- 7
  p <- ggplot() +
    geom_point(data=enrichment_df, aes(x=log10(rep1_enrichment), y=log10(rep2_enrichment)), pch=16, size=1, color='black', alpha=0.5) +
    theme_cowplot(12) +
    theme(plot.title=element_text(size=textsize,face="bold",hjust=0.5),
          axis.title=element_text(size=textsize,face="bold"),
          axis.text=element_text(size=textsize,face="bold"),
          legend.title=element_blank(),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.03, 'in'),
          legend.text=element_text(size=textsize,face="bold"),
          legend.position='top') +
    ggtitle('Italy20 HA') +
    labs(x=expression(bold(log['10']~enrich~'(Replicate 1)')),y=expression(bold(log['10']~enrich~'(Replicate 2)')))
  return (p)
}
count_common_mutation <-function (mutant_list,enrichment_df){
  count_lists <- vector()
  for (i in mutant_list){
    filtered_mutants <- enrichment_df%>%filter(grepl(i, mutants_containing_L194P))
    count_lists <- append(count_lists,nrow(filtered_mutants))}
    count_results <- data.frame(mutant_list,count_lists)
    return(count_results)
}

output1 <- enrichment_correlation(filtered_data, -6)
write.table(output1, file="filtered_L194P_mutant_enrichment.tsv", quote = FALSE, sep='\t', col.names = NA)
p<-graph_correlation(output1)
ggsave(filename='../graph/Replicate_enrich.png',p,height=4,width=3,dpi=600)
l194p <- output1[output1$mutants_containing_L194P == "L194P", ]
improved_enrichment_both_rep <- output1[which(output1$rep1_enrichment>l194p$rep1_enrichment & output1$rep2_enrichment>l194p$rep2_enrichment),]
improved_enrichment_one_rep <- output1[which(output1$rep1_enrichment>l194p$rep1_enrichment | output1$rep2_enrichment>l194p$rep2_enrichment),]
write.table(improved_enrichment_one_rep, file="improved_enrichment_one_rep.tsv", quote = FALSE, sep='\t', col.names = NA)
write.table(improved_enrichment_both_rep, file="improved_enrichment_both_rep.tsv", quote = FALSE, sep='\t', col.names = NA)

count_both_rep<-count_common_mutation(mutant_list, improved_enrichment_both_rep)
count_one_rep<-count_common_mutation(mutant_list, improved_enrichment_one_rep)

write.table(count_both_rep, file="mutation_count_size4.tsv", quote = FALSE, sep='\t', col.names = NA)
write.table(count_one_rep, file="mutation_count_size43.tsv", quote = FALSE, sep='\t', col.names = NA)
