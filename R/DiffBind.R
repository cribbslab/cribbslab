# This script performs the loading of peaksets for ChIP-seq or ATAC-seq using
# the diffbind package. Counts and performed and the output of this is then saved
# as a dba_object.Rdata so it can be loaded in an R session and further analysis
# performed on the data

setwd(".")

suppressPackageStartupMessages(library(optparse))

option_list <- list(
	    make_option(c("--design"), default="must_specify",
	    help="This specifies the design file according to the design file for diffbind."),
	    make_option(c("--output"), default="must_specify",
	    help="This specifies the output Rdata object name.")
	    
)

opt <- parse_args(OptionParser(option_list=option_list))

print("Running with the following options:")
print(opt)


design <- opt$design

suppressPackageStartupMessages(library(DiffBind))

# Load the design file
samples <- read.csv(design)
names(samples)

# Construct the dba object
dba_object <- dba(sampleSheet = samples) 

dba_object

# Plot the correlation heatmap sing peak caller score
plot(dba_object)

#Count reads
dba_object <- dba.count(dba_object, summits=250)

dba_object

# Plot correlation heatmap based on the counts
plot(dba_object)

output <- paste(opt$output, ".Rdata", sep="")
save(dba_object, file=output)
