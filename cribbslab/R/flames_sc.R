#!/usr/bin/env Rscript

# ==============================================================================
# FLAMES single-cell transcript analysis for long-read RNA-seq
#
# Targets the FLAMES 2.x API (Bioconductor >= 3.22 / FLAMES >= 2.4), where the
# pipeline runs end-to-end from FASTQ: barcode demultiplexing -> minimap2 genome
# alignment -> isoform identification -> read realignment -> transcript counting.
# Pre-aligned BAMs are NOT accepted by FLAMES 2.x.
#
# If a barcode whitelist is supplied (e.g. from an upstream BLAZE run), FLAMES is
# configured to demultiplex with flexiplex against that whitelist; otherwise
# FLAMES runs BLAZE internally using the expected cell number.
#
# Usage:
#   Rscript flames_sc.R --fastq <file/dir> --gtf <file> --genome <file> \
#                       --outdir <dir> --threads <int> \
#                       [--barcodes <whitelist>] [--expect_cells <int>]
#
# Output (written by FLAMES into --outdir):
#   - transcript_count.csv.gz, gene_count.csv.gz
#   - isoform_annotated.filtered.gff3, transcript_assembly.fa
#   - align2genome.bam, matched_reads.fastq, ...
# ==============================================================================

suppressPackageStartupMessages({
    library(FLAMES)
    library(optparse)
    library(SingleCellExperiment)
})

# Parse command line arguments
option_list <- list(
    make_option(c("-q", "--fastq"),
                type = "character",
                default = NULL,
                help = "Input FASTQ file (or directory of FASTQs) for one sample",
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
    make_option(c("--barcodes"),
                type = "character",
                default = NULL,
                help = "Optional cell-barcode whitelist (one barcode per line). If given, flexiplex demultiplexing is used.",
                metavar = "FILE"),
    make_option(c("--expect_cells"),
                type = "integer",
                default = 5000,
                help = "Expected number of cells (used when no whitelist is given) [default: %default]",
                metavar = "INT"),
    make_option(c("--min_support"),
                type = "integer",
                default = 3,
                help = "Minimum reads to support a transcript [default: %default]",
                metavar = "INT"),
    make_option(c("--do_discovery"),
                type = "logical",
                default = TRUE,
                help = "Perform novel transcript discovery [default: %default]",
                metavar = "BOOL")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Run FLAMES for single-cell long-read analysis")
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$fastq)) {
    stop("Error: --fastq is required")
}
if (is.null(opt$gtf)) {
    stop("Error: --gtf is required")
}
if (is.null(opt$genome)) {
    stop("Error: --genome is required")
}
if (!file.exists(opt$fastq)) {
    stop(sprintf("Error: --fastq path does not exist: %s", opt$fastq))
}
if (!file.exists(opt$gtf)) {
    stop(sprintf("Error: --gtf does not exist: %s", opt$gtf))
}
if (!file.exists(opt$genome)) {
    stop(sprintf("Error: --genome does not exist: %s", opt$genome))
}

# Create output directory
if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
}

# Decide on the demultiplexing strategy. If a non-empty whitelist is provided,
# demultiplex with flexiplex against it; otherwise run BLAZE with expect_cells.
barcodes_file <- NULL
demultiplexer <- "BLAZE"
if (!is.null(opt$barcodes) &&
        file.exists(opt$barcodes) &&
        file.info(opt$barcodes)$size > 0) {
    barcodes_file <- opt$barcodes
    demultiplexer <- "flexiplex"
}

cat("========================================\n")
cat("FLAMES Single-Cell Analysis\n")
cat("========================================\n")
cat(paste0("FASTQ: ", opt$fastq, "\n"))
cat(paste0("GTF file: ", opt$gtf, "\n"))
cat(paste0("Genome file: ", opt$genome, "\n"))
cat(paste0("Output directory: ", opt$outdir, "\n"))
cat(paste0("Threads: ", opt$threads, "\n"))
cat(paste0("Demultiplexer: ", demultiplexer, "\n"))
if (!is.null(barcodes_file)) {
    cat(paste0("Barcode whitelist: ", barcodes_file, "\n"))
} else {
    cat(paste0("Expected cells: ", opt$expect_cells, "\n"))
}
cat(paste0("Min support reads: ", opt$min_support, "\n"))
cat(paste0("Discovery mode: ", opt$do_discovery, "\n"))
cat("========================================\n\n")

# Build a FLAMES config. create_config() writes a JSON config into outdir and
# returns its path. Parameters are overridden via dot-notation. oarfish
# quantification is disabled to avoid requiring the external 'oarfish' tool.
cat("Configuring FLAMES...\n")
config_file <- create_config(
    opt$outdir,
    type = "sc_3end",
    pipeline_parameters.threads = opt$threads,
    pipeline_parameters.demultiplexer = demultiplexer,
    pipeline_parameters.do_isoform_identification = opt$do_discovery,
    pipeline_parameters.oarfish_quantification = FALSE,
    isoform_parameters.min_sup_cnt = opt$min_support
)
cat(paste0("Configuration written to: ", config_file, "\n\n"))

# Run FLAMES (demultiplex -> align -> isoform ID -> realign -> count)
cat("Running FLAMES analysis...\n")
cat(paste0("Start time: ", Sys.time(), "\n\n"))

sce <- sc_long_pipeline(
    annotation = opt$gtf,
    fastq = opt$fastq,
    genome_fa = opt$genome,
    outdir = opt$outdir,
    barcodes_file = barcodes_file,
    expect_cell_number = if (is.null(barcodes_file)) opt$expect_cells else NULL,
    config_file = config_file
)

cat(paste0("\nEnd time: ", Sys.time(), "\n"))

# sc_long_pipeline returns a SingleCellExperiment on success, or the pipeline
# object if it errored out part-way.
if (!methods::is(sce, "SingleCellExperiment")) {
    stop("FLAMES did not return a SingleCellExperiment - the run failed. ",
         "See the FLAMES messages above and the pipeline.rds in the output ",
         "directory; it can be resumed with FLAMES::resume_FLAMES().")
}

cat("FLAMES analysis complete.\n\n")

# Save the SCE object alongside the FLAMES output files.
cat("Saving output files...\n")
sce_file <- file.path(opt$outdir, "flames_sce.rds")
saveRDS(sce, sce_file)
cat(paste0("SCE object saved to: ", sce_file, "\n"))

# Summary statistics
cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")
cat(paste0("Total transcripts: ", nrow(sce), "\n"))
cat(paste0("Total cells: ", ncol(sce), "\n"))
if (length(assayNames(sce)) > 0) {
    total_counts <- colSums(assay(sce, 1))
    cat("\nCell statistics:\n")
    cat(paste0("  Median counts per cell: ", stats::median(total_counts), "\n"))
    cat(paste0("  Mean counts per cell: ", round(mean(total_counts), 2), "\n"))
}

cat("\n========================================\n")
cat("FLAMES analysis completed successfully!\n")
cat("========================================\n")
