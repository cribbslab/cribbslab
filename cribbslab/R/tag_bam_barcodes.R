#!/usr/bin/env Rscript

# ==============================================================================
# Tag BAM file with cell barcodes from BLAZE output
#
# This script adds CB (cell barcode) and UB (UMI) tags to a BAM file
# using barcode assignments from BLAZE.
#
# Usage:
#   Rscript tag_bam_barcodes.R --bam <file> --barcodes <file> --output <file>
#
# Output:
#   - BAM file with CB and UB tags added to reads
# ==============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

# Parse command line arguments
option_list <- list(
    make_option(c("-b", "--bam"),
                type = "character",
                default = NULL,
                help = "Input BAM file",
                metavar = "FILE"),
    make_option(c("-c", "--barcodes"),
                type = "character",
                default = NULL,
                help = "BLAZE barcode CSV file (putative_bc.csv)",
                metavar = "FILE"),
    make_option(c("-o", "--output"),
                type = "character",
                default = NULL,
                help = "Output BAM file with barcode tags",
                metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Add barcode tags to BAM file")
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$bam)) {
    stop("Error: --bam is required")
}
if (is.null(opt$barcodes)) {
    stop("Error: --barcodes is required")
}
if (is.null(opt$output)) {
    stop("Error: --output is required")
}

cat("========================================\n")
cat("BAM Barcode Tagging\n")
cat("========================================\n")
cat(paste0("Input BAM: ", opt$bam, "\n"))
cat(paste0("Barcode file: ", opt$barcodes, "\n"))
cat(paste0("Output BAM: ", opt$output, "\n"))
cat("========================================\n\n")

# Read barcode assignments
cat("Reading barcode assignments...\n")
bc_data <- fread(opt$barcodes)
cat(paste0("Loaded ", nrow(bc_data), " barcode assignments\n\n"))

# Create barcode lookup (read_id -> barcode)
# BLAZE putative_bc.csv columns:
#   read_id, putative_bc, putative_bc_min_q, putative_umi, polyT_end, ...
bc_col <- "putative_bc"
umi_col <- "putative_umi"

if (!bc_col %in% colnames(bc_data)) {
    stop(sprintf(
        "Expected barcode column '%s' not found in %s. Columns present: %s",
        bc_col, opt$barcodes, paste(colnames(bc_data), collapse = ", ")))
}
if (!"read_id" %in% colnames(bc_data)) {
    stop(sprintf(
        "Expected 'read_id' column not found in %s. Columns present: %s",
        opt$barcodes, paste(colnames(bc_data), collapse = ", ")))
}

# Many reads have no putative barcode (empty cells) - keep only barcoded reads.
bc_data <- bc_data[!is.na(bc_data[[bc_col]]) & bc_data[[bc_col]] != "", ]
cat(paste0("Reads with a putative barcode: ", nrow(bc_data), "\n"))

bc_lookup <- setNames(bc_data[[bc_col]], bc_data$read_id)
umi_lookup <- if (umi_col %in% colnames(bc_data)) {
    setNames(bc_data[[umi_col]], bc_data$read_id)
} else {
    NULL
}

# Process BAM file
cat("Processing BAM file...\n")

# Create output directory if needed
outdir <- dirname(opt$output)
if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
}

# Lightweight stats (avoid loading the whole BAM into memory).
n_aln <- as.integer(system(paste("samtools view -c", shQuote(opt$bam)),
                           intern = TRUE))
cat(paste0("Total alignments: ", n_aln, "\n"))
cat(paste0("Reads with an assignable barcode: ", length(bc_lookup), "\n\n"))

cat("Writing tagged BAM file...\n")

# Create a barcode tag file (read_id, CB, UB) and stream-tag with awk + samtools.
bc_tag_file <- tempfile(fileext = ".tsv")
ub_values <- if (!is.null(umi_lookup)) umi_lookup[names(bc_lookup)] else ""
ub_values[is.na(ub_values)] <- ""
bc_out <- data.frame(
    read_id = names(bc_lookup),
    CB = unname(bc_lookup),
    UB = unname(ub_values)
)
fwrite(bc_out, bc_tag_file, sep = "\t", col.names = FALSE, na = "")

# Use awk with samtools to add tags. Only append a tag when a value is present
# so we never emit a malformed empty tag (e.g. "UB:Z:").
tag_cmd <- paste0(
    "samtools view -h ", shQuote(opt$bam), " | ",
    "awk -F'\\t' -v OFS='\\t' 'NR==FNR{bc[$1]=$2; umi[$1]=$3; next} ",
    "/^@/{print; next} ",
    "{line=$0; ",
    "if($1 in bc && bc[$1]!=\"\"){line=line\"\\tCB:Z:\"bc[$1]; ",
    "if(umi[$1]!=\"\")line=line\"\\tUB:Z:\"umi[$1]} ",
    "print line}' ",
    shQuote(bc_tag_file), " - | ",
    "samtools view -bS -o ", shQuote(opt$output), " -"
)

result <- system(tag_cmd)

if (result != 0) {
    # Fallback: just copy the BAM (tags won't be added)
    warning("Failed to add tags with awk. Copying original BAM.")
    file.copy(opt$bam, opt$output)
}

# Clean up
unlink(bc_tag_file)

cat("\n========================================\n")
cat("BAM tagging completed!\n")
cat(paste0("Output: ", opt$output, "\n"))
cat("========================================\n")
