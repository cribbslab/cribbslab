#!/usr/bin/env Rscript

# ==============================================================================
# Tag BAM file with cell barcodes from BLAZE output
#
# Applies BLAZE step-3 read-to-whitelist assignment (error-corrected CB/UB)
# before writing CB and UB tags. Using raw putative_bc values causes many
# reads to miss subset_cells when CB does not exactly match the whitelist.
#
# Usage:
#   Rscript tag_bam_barcodes.R --bam <file> --barcodes <putative_bc.csv> \
#       [--whitelist <whitelist.csv>] --output <file>
#
# Output:
#   - BAM file with CB and UB tags added to reads
# ==============================================================================

suppressPackageStartupMessages({
    library(optparse)
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
                help = "BLAZE putative_bc.csv",
                metavar = "FILE"),
    make_option(c("-w", "--whitelist"),
                type = "character",
                default = NULL,
                help = paste0(
                    "BLAZE whitelist.csv for corrected assignments ",
                    "(default: infer from --barcodes path)"
                ),
                metavar = "FILE"),
    make_option(c("-o", "--output"),
                type = "character",
                default = NULL,
                help = "Output BAM file with barcode tags",
                metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Add corrected BLAZE barcode tags to BAM file")
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
if (!file.exists(opt$barcodes)) {
    stop(sprintf("Error: barcode file not found: %s", opt$barcodes))
}

# Resolve sibling whitelist when not supplied explicitly.
whitelist_path <- opt$whitelist
if (is.null(whitelist_path) || !nzchar(whitelist_path)) {
    whitelist_path <- sub("_putative_bc\\.csv$", "_whitelist.csv", opt$barcodes)
}

cat("========================================\n")
cat("BAM Barcode Tagging (BLAZE-corrected)\n")
cat("========================================\n")
cat(paste0("Input BAM: ", opt$bam, "\n"))
cat(paste0("Putative BC: ", opt$barcodes, "\n"))
cat(paste0("Whitelist: ", whitelist_path, "\n"))
cat(paste0("Output BAM: ", opt$output, "\n"))
cat("========================================\n\n")

if (!file.exists(whitelist_path)) {
    stop(sprintf(
        paste0(
            "Error: BLAZE whitelist not found: %s\n",
            "Pass --whitelist or ensure BLAZE completed step 2."
        ),
        whitelist_path
    ))
}

# Locate blaze_assign_barcodes.py relative to this script (../python/).
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg)) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
    script_dir <- getwd()
}
assign_script <- normalizePath(
    file.path(script_dir, "..", "python", "blaze_assign_barcodes.py"),
    mustWork = TRUE
)

python_bin <- Sys.which("python3")
if (!nzchar(python_bin)) {
    python_bin <- Sys.which("python")
}
if (!nzchar(python_bin)) {
    stop("Error: python3/python not found on PATH")
}

assignments_tsv <- tempfile(fileext = ".blaze_assignments.tsv")
assign_cmd <- paste(
    shQuote(python_bin),
    shQuote(assign_script),
    "--putative-bc", shQuote(opt$barcodes),
    "--whitelist", shQuote(whitelist_path),
    "--output", shQuote(assignments_tsv)
)
cat("Running BLAZE read assignment...\n")
cat(paste0("Command: ", assign_cmd, "\n\n"))
assign_status <- system(assign_cmd)
if (assign_status != 0) {
    stop("BLAZE read assignment failed")
}

n_assigned <- length(readLines(assignments_tsv))
cat(paste0("Reads with corrected CB assignment: ", n_assigned, "\n\n"))

# Process BAM file
cat("Processing BAM file...\n")

outdir <- dirname(opt$output)
if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
}

n_aln <- as.integer(system(paste("samtools view -c", shQuote(opt$bam)),
                           intern = TRUE))
cat(paste0("Total alignments: ", n_aln, "\n\n"))

cat("Writing tagged BAM file...\n")

# assignments_tsv columns: read_id, CB, UB
tag_cmd <- paste0(
    "samtools view -h ", shQuote(opt$bam), " | ",
    "awk -F'\\t' -v OFS='\\t' 'NR==FNR{bc[$1]=$2; umi[$1]=$3; next} ",
    "/^@/{print; next} ",
    "{line=$0; ",
    "if($1 in bc && bc[$1]!=\"\"){line=line\"\\tCB:Z:\"bc[$1]; ",
    "if(umi[$1]!=\"\")line=line\"\\tUB:Z:\"umi[$1]} ",
    "print line}' ",
    shQuote(assignments_tsv), " - | ",
    "samtools view -bS -o ", shQuote(opt$output), " -"
)

result <- system(tag_cmd)

if (result != 0) {
    warning("Failed to add tags with awk. Copying original BAM.")
    file.copy(opt$bam, opt$output)
}

unlink(assignments_tsv)

cat("\n========================================\n")
cat("BAM tagging completed!\n")
cat(paste0("Output: ", opt$output, "\n"))
cat("========================================\n")
