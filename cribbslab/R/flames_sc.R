#!/usr/bin/env Rscript

# ==============================================================================
# FLAMES single-cell transcript analysis for long-read RNA-seq
#
# This script runs FLAMES on single-cell long-read data to:
# 1. Discover novel transcripts/isoforms
# 2. Quantify transcript expression per cell
# 3. Generate gene-level and transcript-level count matrices
# 4. Generate spliced/unspliced matrices for RNA velocity
#
# Usage:
#   Rscript flames_sc.R --bam_list <file> --gtf <file> --genome <file> \
#                       --outdir <dir> --threads <int>
#
# Output:
#   - transcript_count.csv.gz: Transcript-level counts per cell
#   - gene_count.csv.gz: Gene-level counts per cell
#   - isoform_annotated.filtered.gff3: Novel isoform annotations
#   - transcript_assembly.fa: Transcript sequences
# ==============================================================================

suppressPackageStartupMessages({
    library(FLAMES)
    library(optparse)
    library(SingleCellExperiment)
})

# Parse command line arguments
option_list <- list(
    make_option(c("-b", "--bam_list"),
                type = "character",
                default = NULL,
                help = "File containing list of BAM file paths (one per line)",
                metavar = "FILE"),
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
                default = "flames_output",
                help = "Output directory [default: %default]",
                metavar = "DIR"),
    make_option(c("-t", "--threads"),
                type = "integer",
                default = 8,
                help = "Number of threads [default: %default]",
                metavar = "INT"),
    make_option(c("--min_support"),
                type = "integer",
                default = 3,
                help = "Minimum reads to support a transcript [default: %default]",
                metavar = "INT"),
    make_option(c("--min_cells"),
                type = "integer",
                default = 3,
                help = "Minimum cells expressing a transcript [default: %default]",
                metavar = "INT"),
    make_option(c("--do_discovery"),
                type = "logical",
                default = TRUE,
                action = "store_true",
                help = "Perform novel transcript discovery [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Run FLAMES for single-cell long-read analysis")
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$bam_list)) {
    stop("Error: --bam_list is required")
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
cat("FLAMES Single-Cell Analysis\n")
cat("========================================\n")
cat(paste0("BAM list: ", opt$bam_list, "\n"))
cat(paste0("GTF file: ", opt$gtf, "\n"))
cat(paste0("Genome file: ", opt$genome, "\n"))
cat(paste0("Output directory: ", opt$outdir, "\n"))
cat(paste0("Threads: ", opt$threads, "\n"))
cat(paste0("Min support reads: ", opt$min_support, "\n"))
cat(paste0("Min cells: ", opt$min_cells, "\n"))
cat(paste0("Discovery mode: ", opt$do_discovery, "\n"))
cat("========================================\n\n")

# Read BAM file list
bam_files <- readLines(opt$bam_list)
bam_files <- bam_files[bam_files != ""]

cat(paste0("Found ", length(bam_files), " BAM files:\n"))
for (f in bam_files) {
    cat(paste0("  - ", basename(f), "\n"))
}
cat("\n")

# Configure FLAMES
cat("Configuring FLAMES...\n")

config <- jsonlite::fromJSON(
    system.file("extdata", "config_sclr_nanopore_default.json", package = "FLAMES")
)

# Update config with user parameters
config$pipeline_parameters$do_isoform_identification <- opt$do_discovery
config$isoform_parameters$min_sup_cnt <- opt$min_support
config$alignment_parameters$threads <- opt$threads

# Save config
config_file <- file.path(opt$outdir, "flames_config.json")
jsonlite::write_json(config, config_file, auto_unbox = TRUE, pretty = TRUE)

cat("Configuration saved to: ", config_file, "\n\n")

# Run FLAMES
cat("Running FLAMES analysis...\n")
cat(paste0("Start time: ", Sys.time(), "\n\n"))

# Create FLAMES SingleCellExperiment
sce <- sc_long_pipeline(
    annotation = opt$gtf,
    fastq = NULL,  # We're using pre-aligned BAMs
    genome_fa = opt$genome,
    outdir = opt$outdir,
    bam_path = bam_files,
    config_file = config_file
)

cat(paste0("\nEnd time: ", Sys.time(), "\n"))
cat("FLAMES analysis complete.\n\n")

# Save outputs
cat("Saving output files...\n")

# Save SCE object
sce_file <- file.path(opt$outdir, "flames_sce.rds")
saveRDS(sce, sce_file)
cat(paste0("SCE object saved to: ", sce_file, "\n"))

# Extract and save count matrices
if (!is.null(sce)) {
    # Transcript counts
    transcript_counts <- counts(sce)
    transcript_file <- file.path(opt$outdir, "transcript_counts.csv.gz")
    write.csv(as.matrix(transcript_counts), gzfile(transcript_file))
    cat(paste0("Transcript counts saved to: ", transcript_file, "\n"))
    
    # Gene counts (aggregate transcripts to genes)
    if ("gene_id" %in% colnames(rowData(sce))) {
        gene_counts <- rowsum(as.matrix(transcript_counts), rowData(sce)$gene_id)
        gene_file <- file.path(opt$outdir, "gene_counts.csv.gz")
        write.csv(gene_counts, gzfile(gene_file))
        cat(paste0("Gene counts saved to: ", gene_file, "\n"))
    }
}

# Summary statistics
cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")

if (!is.null(sce)) {
    cat(paste0("Total transcripts: ", nrow(sce), "\n"))
    cat(paste0("Total cells: ", ncol(sce), "\n"))
    
    # Transcript class summary
    if ("transcript_class" %in% colnames(rowData(sce))) {
        class_counts <- table(rowData(sce)$transcript_class)
        cat("\nTranscript classes:\n")
        for (class_name in names(class_counts)) {
            cat(paste0("  ", class_name, ": ", class_counts[class_name], "\n"))
        }
    }
    
    # Cell statistics
    total_counts <- colSums(counts(sce))
    cat("\nCell statistics:\n")
    cat(paste0("  Median counts per cell: ", median(total_counts), "\n"))
    cat(paste0("  Mean counts per cell: ", round(mean(total_counts), 2), "\n"))
    cat(paste0("  Min counts per cell: ", min(total_counts), "\n"))
    cat(paste0("  Max counts per cell: ", max(total_counts), "\n"))
}

cat("\n========================================\n")
cat("FLAMES analysis completed successfully!\n")
cat("========================================\n")

# Save session info
session_file <- file.path(opt$outdir, "flames_session_info.txt")
writeLines(capture.output(sessionInfo()), session_file)
cat(paste0("\nSession info saved to: ", session_file, "\n"))
