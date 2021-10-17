# This script parses the output from MACS3 and then annotates the regions, it will generate 

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ChIPpeakAnno))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

option_list <- list(
	    make_option(c("--xls"), default="must_specify",
	    help="This specifys the xls output from MACS3"),
	    make_option(c("--name"), default="must_specify",
	    help="This specifys the name of the file "))

opt <- parse_args(OptionParser(option_list=option_list))

print("Running with the following options:")
print(opt)


macsOutput <- ChIPpeakAnno::toGRanges(opt$xls, format="MACS2")


suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

ucsc.hg38.knownGene <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
macs.anno <- annotatePeakInBatch(macsOutput, 
                                     AnnotationData=ucsc.hg38.knownGene)

print(macs.anno)
macs.anno <- addGeneIDs(annotatedPeak=macs.anno, 
                            orgAnn="org.Hs.eg.db", 
                            feature_id_type="entrez_id",
                            IDs2Add="symbol")

print(macs.anno)



df <- data.frame(seqnames=seqnames(macs.anno),
  starts=start(macs.anno)-1,
  ends=end(macs.anno),
  names=c(rep(".", length(macs.anno))),
  scores=c(rep(".", length(macs.anno))),
  strands=strand(macs.anno),
  symbol=macs.anno$symbol)

bedfile = paste0(opt$name, "_annotatePeaks.bed")

write.table(df, file=bedfile, quote=F, sep="\t", row.names=F, col.names=F)


filename = paste0(opt$name,"_ChIPpeakAnno.Rdata")

save(macs.anno, file=filename)


