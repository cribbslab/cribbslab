# Coverage plot test script

library("tidyverse")
library("ggpubr")

# Generalised version, need to make standardised 
setwd("/ifs/research-groups/botnar/proj025/analyses/test_trna_3")

# Just using one sample, would need to make multiple plots for multiple files (dynamically scale)
pileup <- read_tsv(file = "post_mapping_bams.dir/test_pileup.tsv", skip = 3, col_names = TRUE)

# Only subset of 5 clusters. Find out what top 50 means, top 50 differentially expressed--> where is this info located
small_pileup <- pileup %>% dplyr::filter(`#CHROM` %in% head(unique(pileup$`#CHROM`),60))

subset_pileup <- dplyr::select(small_pileup, c(`#CHROM`, POS, REF, ALT)) %>% dplyr::rename(cluster = `#CHROM`, Base = REF)
subset_pileup$Base <- as.factor(subset_pileup$Base)
subset_pileup$ALT <- as.factor(subset_pileup$ALT)

coverage_plot <- ggplot(data = subset_pileup, aes(x= POS, y = cluster, color = Base)) + 
  geom_point() +theme_bw() + labs(x= "Position in tRNA", y = "Cluster")
  
set_palette(coverage_plot, "lancet")

# Need to convert to AA / anticodon / 5' vs 3' identifiers?? e.g. Gly-GCC-5'