library(ATACseqQC)
library(BSgenome.Hsapiens.UCSC.hg19)

bamfile <- "/Users/adamc/Documents/cgat-developers/cgat-apps/tests/bam2bam.py/paired.bam"
bamfile.labels <- gsub(".bam", "", basename(bamfile))

source(system.file("extdata", "IGVSnapshot.R", package = "ATACseqQC"))

estimateLibComplexity(readsDupFreq(bamfile))

fragSize <- fragSizeDist(bamfile, bamfile.labels)


## bamfile tags to be read in
tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
## files will be output into outPath
outPath <- "splited"
dir.create(outPath)

seqlev <- "chr1" ## subsample data for quick run
which <- as(seqinfo(Hsapiens)[seqlev], "GRanges")
gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
gal1 <- shiftGAlignmentsList(gal)
shiftedBamfile <- file.path(outPath, "shifted.bam")
export(gal1, shiftedBamfile)