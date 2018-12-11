#!/usr/bin/Rscript

library("ggplot2")
library("tidyverse")
library("optparse")
library("ggpubr")

option_list <- list(
	    make_option(c("--input"), default="must_specify",
	    help="this specifiys the input output from count_features"),
	    make_option(c("--output"), default="must_specify",
	    help="this specifiys the input output from count_features")
)
#OptionParser(option_list=option_list)
opt <- parse_args(OptionParser(option_list=option_list))

print("Running with following options:")
print(opt)

#test_in <- "featurecounts.dir/test/test.feature_small.tsv:featurecounts.dir/test_2/test_2.feature_small.tsv"
#test_out <- "plots.dir/feature_count.png"

#input <- unlist(strsplit(test_in, ":"))
input <- unlist(strsplit(opt$input, ":"))
output <- opt$output

num_samples = length(input)
feature_count_tibble_all_samples = NULL

# Iterate over each sample
for(i in 1:num_samples){
  path <- input[i]
  feature_count_tibble <-read_tsv(path, col_names = FALSE, skip = 2)
  colnames(feature_count_tibble) <- c("GeneID", "NumberMapped")
  # Get sample name 
  sample <- stringr::str_match(path, "/(.*?)/")[2]
  feature_count_tibble$Sample <- sample
  total_mapped <- sum(feature_count_tibble$NumberMapped)
  feature_count_tibble$Percentages <- 100*(feature_count_tibble$NumberMapped / total_mapped)
  feature_count_tibble_all_samples = bind_rows(feature_count_tibble_all_samples,feature_count_tibble)
}
fc_tibble_all_samples_nonzero <- dplyr::filter(feature_count_tibble_all_samples, NumberMapped != 0)

feature_count_plot <- ggplot(data=fc_tibble_all_samples_nonzero, aes(x=Sample, y=Percentages, fill=`GeneID`)) +
  geom_bar(stat="identity") + coord_flip() 
  
feature_count_plot <- set_palette(feature_count_plot, "lancet")

#print(feature_count_plot)

ggsave(filename = output, plot = feature_count_plot)


file_paths = list.files(pattern='_small.tsv$', recursive=TRUE, full.names = TRUE)
num_samples <- length(file_paths)
feature_count_tibble_all_samples = NULL

for(i in 1:num_samples){
  path <- file_paths[i]
  feature_count_tibble <-read_tsv(path, col_names = FALSE, skip = 2)
  colnames(feature_count_tibble) <- c("GeneID", "NumberMapped")
  # Get sample name 
  sample <- stringr::str_match(path, "/(.*?)/")[2]
  feature_count_tibble$Sample <- sample
  total_mapped <- sum(feature_count_tibble$NumberMapped)
  feature_count_tibble$Percentages <- 100*(feature_count_tibble$NumberMapped / total_mapped)
  feature_count_tibble_all_samples = bind_rows(feature_count_tibble_all_samples,feature_count_tibble)
}
fc_tibble_all_samples_nonzero <- dplyr::filter(feature_count_tibble_all_samples, NumberMapped != 0)

feature_count_plot <- ggplot(data=fc_tibble_all_samples_nonzero, aes(x=Sample, y=Percentages, fill=`GeneID`)) +
  geom_bar(stat="identity") + coord_flip() 

feature_count_plot <- set_palette(feature_count_plot, "lancet")

print(feature_count_plot)
