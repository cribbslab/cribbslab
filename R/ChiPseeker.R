# This script parses the narrowPeaks from control and treatment file output from macs2 and saves
# the data as two seperate R objects.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
	    make_option(c("--control"), default="must_specify",
	    help="This specifys the directory and wildcard that specifies where the control narrowPeaks files are located"),
	    make_option(c("--treatment"), default="must_specify",
	    help="This specifys the diredctory and wildcard that specifies where the treatment narrowPeaks files are located")
)

opt <- parse_args(OptionParser(option_list=option_list))

print("Running with the following options:")
print(opt)


dir_control <- opt$control
dir_treatment <- opt$treatment

suppressPackageStartupMessages(library(rtracklayer))

extracols <- c(signalValue = "numeric", pvalue = "numeric", qValue = "numeric", peak = "integer")

peaks_control <- lapply(Sys.glob(dir_control), function(i){
  import(i, format= "BED", extraCols= extracols)})

peaks_treatment <- lapply(Sys.glob(dir_treatment), function(i){
  import(i, format= "BED", extraCols= extracols)})

save(peaks_control, peaks_treatment, file="ChIPseeker.Rdata")


