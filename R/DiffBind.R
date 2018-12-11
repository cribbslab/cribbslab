# This script performs the loading of peaksets from the 

setwd(".")

suppressPackageStartupMessages(library(optparse))

option_list <- list(
	    make_option(c("--design"), default="must_specify",
	    help="This specifies the design file according to the design file for diffbind.")
)

opt <- parse_args(OptionParser(option_list=option_list))

print("Running with the following options:")
print(opt)


design <- opt$design


save(peaks_control, peaks_treatment, file="ChIPseeker.Rdata")


