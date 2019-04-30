# Coverage plot test script

library("tidyverse")
library("ggpubr")

# Generalised version, need to make standardised 
setwd("/ifs/research-groups/botnar/proj025/analyses/test_april/")

# Just using one sample, would need to make multiple plots for multiple files (dynamically scale)
pileup <- read_tsv(file = "post_mapping_bams.dir/test_pileup.tsv", skip = 3, col_names = TRUE)

# Only subset of 5 clusters. Find out what top 50 means, top 50 differentially expressed--> where is this info located
small_pileup <- pileup %>% dplyr::filter(`#CHROM` %in% head(unique(pileup$`#CHROM`),60))

subset_pileup <- dplyr::select(small_pileup, c(`#CHROM`, POS, REF, ALT)) %>% dplyr::rename(cluster = `#CHROM`, Base = REF)
subset_pileup$Base <- as.factor(subset_pileup$Base)
subset_pileup$ALT <- as.factor(subset_pileup$ALT)

# Rank clusters by abundance
idx_stats <- read_tsv(file = "post_mapping_bams.dir/test_trna.idxstats", col_names = FALSE)
colnames(idx_stats) <- c("cluster", "length", "abundance", "not_sure")

idx_stats_filt <-  arrange(idx_stats, desc(abundance)) %>% top_n(50, abundance)
top_clusters <- idx_stats_filt$cluster

subset_pileup_top <- subset_pileup %>% filter(cluster %in% top_clusters) 
subset_pileup_top$cluster <- factor(subset_pileup_top$cluster, levels = rev(top_clusters))

set_palette(coverage_plot, "lancet")
coverage_plot <- ggplot(data = subset_pileup_top, aes(x= POS, y = cluster, color = Base)) + 
  geom_point() +theme_bw() + labs(x= "Position in tRNA", y = "Cluster")
print(coverage_plot)
