#!/usr/bin/env Rscript

# ==============================================================================
# Bambu transcript discovery and quantification for long-read RNA-seq
#
# This script runs bambu on aligned BAM files to:
# 1. Discover novel transcripts
# 2. Quantify transcript expression
#
# Usage:
#   Rscript bambu.R --bam_dir <dir> --gtf <file> --genome <file> --outdir <dir> \
#                   --prefix <str> --threads <int>
#
# Output:
#   - counts_gene.txt: Gene-level counts
#   - counts_transcript.txt: Transcript-level counts
#   - extended_annotations.gtf: GTF with novel transcripts
#   - CPM.txt: Counts per million
#   - fullLengthCounts.txt: Full-length transcript counts
# ==============================================================================

suppressPackageStartupMessages({
    library(bambu)
    library(optparse)
})

# Parse command line arguments
option_list <- list(
    make_option(c("-b", "--bam_dir"),
                type = "character",
                default = NULL,
                help = "Directory containing BAM files",
                metavar = "DIR"),
    make_option(c("-g", "--gtf"),
                type = "character",
                default = NULL,
                help = "Path to annotation GTF file",
                metavar = "FILE"),
    make_option(c("-f", "--genome"),
                type = "character",
                default = NULL,
                help = "Path to reference genome FASTA file",
                metavar = "FILE"),
    make_option(c("-o", "--outdir"),
                type = "character",
                default = "bambu_output",
                help = "Output directory [default: %default]",
                metavar = "DIR"),
    make_option(c("-p", "--prefix"),
                type = "character",
                default = "sample",
                help = "Output file prefix [default: %default]",
                metavar = "STR"),
    make_option(c("-t", "--threads"),
                type = "integer",
                default = 8,
                help = "Number of threads [default: %default]",
                metavar = "INT"),
    make_option(c("-a", "--annotations_rds"),
                type = "character",
                default = NULL,
                help = "Optional: Path to pre-prepared annotations RDS file",
                metavar = "FILE"),
    make_option(c("--discovery"),
                type = "logical",
                default = TRUE,
                action = "store_true",
                help = "Enable novel transcript discovery [default: %default]"),
    make_option(c("--no_discovery"),
                action = "store_false",
                dest = "discovery",
                help = "Disable novel transcript discovery"),
    make_option(c("--save_annotations"),
                type = "logical",
                default = TRUE,
                action = "store_true",
                help = "Save prepared annotations as RDS for future use [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Run bambu for transcript discovery and quantification")
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$bam_dir)) {
    stop("Error: --bam_dir is required")
}
if (is.null(opt$gtf)) {
    stop("Error: --gtf is required")
}
if (is.null(opt$genome)) {
    stop("Error: --genome is required")
}

# Create output directory
if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
}

cat("========================================\n")
cat("Bambu Transcript Analysis\n")
cat("========================================\n")
cat(paste0("BAM directory: ", opt$bam_dir, "\n"))
cat(paste0("GTF file: ", opt$gtf, "\n"))
cat(paste0("Genome file: ", opt$genome, "\n"))
cat(paste0("Output directory: ", opt$outdir, "\n"))
cat(paste0("Prefix: ", opt$prefix, "\n"))
cat(paste0("Threads: ", opt$threads, "\n"))
cat(paste0("Discovery mode: ", opt$discovery, "\n"))
cat("========================================\n\n")

# Get BAM files
bam_files <- list.files(opt$bam_dir, pattern = "\\.bam$", full.names = TRUE)

if (length(bam_files) == 0) {
    stop(paste0("Error: No BAM files found in ", opt$bam_dir))
}

cat(paste0("Found ", length(bam_files), " BAM files:\n"))
for (f in bam_files) {
    cat(paste0("  - ", basename(f), "\n"))
}
cat("\n")

# Prepare or load annotations
if (!is.null(opt$annotations_rds) && file.exists(opt$annotations_rds)) {
    cat("Loading pre-prepared annotations from RDS file...\n")
    annotations <- readRDS(opt$annotations_rds)
} else {
    cat("Preparing annotations from GTF file...\n")
    annotations <- prepareAnnotations(opt$gtf)
    
    # Save annotations for future use
    if (opt$save_annotations) {
        annotations_rds <- file.path(opt$outdir, paste0(opt$prefix, "_annotations.rds"))
        cat(paste0("Saving annotations to: ", annotations_rds, "\n"))
        saveRDS(annotations, annotations_rds)
    }
}

cat("\nAnnotations prepared successfully.\n\n")

# Run bambu
cat("Running bambu...\n")
cat(paste0("Start time: ", Sys.time(), "\n\n"))

bambu_result <- bambu(
    reads = bam_files,
    annotations = annotations,
    genome = opt$genome,
    discovery = opt$discovery,
    verbose = TRUE,
    ncore = opt$threads
)

cat(paste0("\nEnd time: ", Sys.time(), "\n"))
cat("Bambu analysis complete.\n\n")

# Write output
cat("Writing output files...\n")

writeBambuOutput(
    bambu_result,
    path = opt$outdir,
    prefix = paste0(opt$prefix, "_")
)

cat("\nOutput files written to: ", opt$outdir, "\n")
cat("Files created:\n")
output_files <- list.files(opt$outdir, pattern = opt$prefix)
for (f in output_files) {
    cat(paste0("  - ", f, "\n"))
}

# Summary statistics
cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")

# Get counts
se <- bambu_result
counts <- assays(se)$counts

cat(paste0("Total genes: ", nrow(counts), "\n"))
cat(paste0("Total samples: ", ncol(counts), "\n"))

# Calculate some basic stats
total_counts <- colSums(counts)
cat("\nTotal counts per sample:\n")
for (i in seq_along(total_counts)) {
    cat(paste0("  ", names(total_counts)[i], ": ", format(total_counts[i], big.mark = ","), "\n"))
}

# Novel transcript discovery summary (if enabled)
if (opt$discovery) {
    rowdata <- rowData(se)
    if ("newTxClass" %in% colnames(rowdata)) {
        novel_counts <- table(rowdata$newTxClass)
        cat("\nTranscript classes:\n")
        for (class_name in names(novel_counts)) {
            cat(paste0("  ", class_name, ": ", novel_counts[class_name], "\n"))
        }
    }
}

cat("\n========================================\n")
cat("Bambu analysis completed successfully!\n")
cat("========================================\n")

# Save session info for reproducibility
session_file <- file.path(opt$outdir, paste0(opt$prefix, "_session_info.txt"))
writeLines(capture.output(sessionInfo()), session_file)
cat(paste0("\nSession info saved to: ", session_file, "\n"))
