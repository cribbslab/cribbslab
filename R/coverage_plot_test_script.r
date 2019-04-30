# Coverage plot test script

library("tidyverse")
library("ggpubr")
library("optparse")


option_list <- list(
  make_option(c("--input"), default="must_specify",
              help="*.pileup.tsv file. Output from create_coverage"),
  make_option(c("--output"), default="must_specify",
              help="file for plot")
)
#OptionParser(option_list=option_list)
opt <- parse_args(OptionParser(option_list=option_list))

output_file_path <- opt$output
pileup_file_path <- opt$input

sample <- gsub("_pileup.tsv", "", basename(pileup_file_path))

# To test
# setwd("/ifs/research-groups/botnar/proj025/analyses/test_april/")

# Just using one sample, would need to make multiple plots for multiple files (dynamically scale)
pileup <- read_tsv(file = pileup_file_path, skip = 3, col_names = TRUE)

# Only subset of 5 clusters. Find out what top 50 means, top 50 differentially expressed--> where is this info located
small_pileup <- pileup %>% dplyr::filter(`#CHROM` %in% head(unique(pileup$`#CHROM`),60))

subset_pileup <- dplyr::select(small_pileup, c(`#CHROM`, POS, REF, ALT)) %>% dplyr::rename(cluster = `#CHROM`, Base = REF)
subset_pileup$Base <- as.factor(subset_pileup$Base)
subset_pileup$ALT <- as.factor(subset_pileup$ALT)

# Rank clusters by abundance
idx_stats <- read_tsv(file = paste("post_mapping_bams.dir/", sample, "_trna.idxstats", sep = ""), col_names = FALSE)
colnames(idx_stats) <- c("cluster", "length", "abundance", "not_sure")

idx_stats_filt <-  arrange(idx_stats, desc(abundance)) %>% top_n(50, abundance)
top_clusters <- idx_stats_filt$cluster

subset_pileup_top <- subset_pileup %>% filter(cluster %in% top_clusters) 

## 
# Abundance per base
table_full <- tibble(cluster = factor(), POS = double(), percent = double())
for(cluster in top_clusters){
  file_name = sub('-', '.csv', cluster)
  path = paste("tRNA-end-site.dir/", sample, ".dir/", file_name, sep = "")
  table <- read_csv(path, col_names = c("POS", "eh", "percent"), skip =1)
  table$cluster <- cluster
  table <- table %>% select(cluster, POS, percent)
  table_full <- rbind(table_full, table)
}

lj <- left_join(table_full, subset_pileup_top , c("cluster" = "cluster", "POS" = "POS") )
lj$cluster <- factor(lj$cluster, levels = rev(top_clusters))

#coverage_plot <- ggplot(data = lj, aes(x= POS, y = cluster, color = Base, alpha = percent)) + 
 # geom_point() + theme_bw() + labs(x= "Position in tRNA", y = "Cluster")

#print(coverage_plot)
png(output_file_path)
ggplot(data = lj, aes(x= POS, y = cluster, color = Base, alpha = percent)) + 
  geom_point() + theme_bw() + labs(x= "Position in tRNA", y = "Cluster")

dev.off()

#rj <- right_join(table_full, subset_pileup_top , c("cluster" = "cluster", "POS" = "POS") )
