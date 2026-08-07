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
#       [--whitelist <whitelist.csv>] --output <file> [--threads <int>]
#
# Output:
#   - BAM file with CB and UB tags added to reads
#   - Sibling assignments TSV (reused on rerun) and progress log
# ==============================================================================

suppressPackageStartupMessages({
    library(optparse)
})

flush_msg <- function(...) {
    cat(...)
    flush.console()
}

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
                metavar = "FILE"),
    make_option(c("-t", "--threads"),
                type = "integer",
                default = 4,
                help = "Threads for BLAZE read assignment [default: %default]",
                metavar = "INT")
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
if (!file.exists(opt$bam)) {
    stop(sprintf("Error: BAM file not found: %s", opt$bam))
}

# Resolve sibling whitelist when not supplied explicitly.
whitelist_path <- opt$whitelist
if (is.null(whitelist_path) || !nzchar(whitelist_path)) {
    whitelist_path <- sub("_putative_bc\\.csv$", "_whitelist.csv", opt$barcodes)
}

flush_msg("========================================\n")
flush_msg("BAM Barcode Tagging (BLAZE-corrected)\n")
flush_msg("========================================\n")
flush_msg(paste0("Input BAM: ", opt$bam, "\n"))
flush_msg(paste0("Putative BC: ", opt$barcodes, "\n"))
flush_msg(paste0("Whitelist: ", whitelist_path, "\n"))
flush_msg(paste0("Output BAM: ", opt$output, "\n"))
flush_msg(paste0("Assignment threads: ", opt$threads, "\n"))
flush_msg("========================================\n\n")
flush_msg(paste0(
    "Note: cgatcore may hide this log until the job finishes.\n",
    "Watch the progress log below with: tail -f <progress.log>\n\n"
))

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

outdir <- dirname(opt$output)
if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
}

# Persist assignments next to the BAM so interrupted jobs can resume without
# redoing the expensive ED matching step.
sample_stem <- sub("\\.tagged\\.bam$", "", basename(opt$output))
assignments_tsv <- file.path(outdir, paste0(sample_stem, ".blaze_assignments.tsv"))
progress_log <- file.path(outdir, paste0(sample_stem, ".tag_progress.log"))

if (file.exists(assignments_tsv) && file.info(assignments_tsv)$size > 0) {
    flush_msg(paste0(
        "Reusing existing assignments: ", assignments_tsv, "\n",
        "(delete this file to force re-assignment)\n\n"
    ))
} else {
    assign_cmd <- paste(
        shQuote(python_bin),
        shQuote(assign_script),
        "--putative-bc", shQuote(opt$barcodes),
        "--whitelist", shQuote(whitelist_path),
        "--output", shQuote(assignments_tsv),
        "--threads", as.integer(opt$threads),
        "--progress-log", shQuote(progress_log)
    )
    flush_msg("Running BLAZE read assignment (this is the slow step)...\n")
    flush_msg(paste0("Command: ", assign_cmd, "\n"))
    flush_msg(paste0("Progress log: ", progress_log, "\n\n"))
    assign_status <- system(assign_cmd)
    if (assign_status != 0) {
        stop("BLAZE read assignment failed")
    }
}

n_assigned <- as.integer(system(
    paste("wc -l <", shQuote(assignments_tsv)),
    intern = TRUE
))
flush_msg(paste0("Reads with corrected CB assignment: ", n_assigned, "\n\n"))

flush_msg("Writing tagged BAM (streaming; no pre-count of alignments)...\n")
flush_msg(paste0(
    "This can take a long time on large BAMs with no intermediate output.\n",
    "Watch output file size: ls -lh ", shQuote(opt$output), "\n\n"
))

# assignments_tsv columns: read_id, CB, UB
# Load assignments into awk memory, then stream the BAM once.
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
    stop("Failed to write tagged BAM")
}

if (!file.exists(opt$output) || file.info(opt$output)$size == 0) {
    stop(sprintf("Tagged BAM missing or empty: %s", opt$output))
}

flush_msg("\n========================================\n")
flush_msg("BAM tagging completed!\n")
flush_msg(paste0("Output: ", opt$output, "\n"))
flush_msg(paste0("Assignments: ", assignments_tsv, "\n"))
flush_msg("========================================\n")
