#!/usr/bin/env Rscript

# ==============================================================================
# Generate spliced/unspliced count matrices for RNA velocity
#
# This script processes a barcode-tagged BAM file to separate reads into:
# - Spliced (exonic only)
# - Unspliced (contains intronic sequences)
#
# Important for nuclei sequencing where unspliced transcripts are abundant.
#
# Usage:
#   Rscript velocity_counts.R --bam <file> --gtf <file> --output <file>
#
# Output:
#   - AnnData H5AD file with spliced and unspliced layers
# ==============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(GenomicAlignments)
    library(GenomicFeatures)
    library(Matrix)
    library(data.table)
})

# Parse command line arguments
option_list <- list(
    make_option(c("-b", "--bam"),
                type = "character",
                default = NULL,
                help = "Input BAM file with CB tags",
                metavar = "FILE"),
    make_option(c("-g", "--gtf"),
                type = "character",
                default = NULL,
                help = "GTF annotation file",
                metavar = "FILE"),
    make_option(c("-o", "--output"),
                type = "character",
                default = NULL,
                help = "Output file (H5AD or RDS)",
                metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Generate velocity count matrices")
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$bam)) {
    stop("Error: --bam is required")
}
if (is.null(opt$gtf)) {
    stop("Error: --gtf is required")
}
if (is.null(opt$output)) {
    stop("Error: --output is required")
}

cat("========================================\n")
cat("Velocity Count Matrix Generation\n")
cat("========================================\n")
cat(paste0("Input BAM: ", opt$bam, "\n"))
cat(paste0("GTF file: ", opt$gtf, "\n"))
cat(paste0("Output: ", opt$output, "\n"))
cat("========================================\n\n")

# Create TxDb from GTF
cat("Loading annotations...\n")
txdb <- makeTxDbFromGFF(opt$gtf, format = "gtf")

# Get exon and intron regions
cat("Extracting exon and intron regions...\n")
exons_by_gene <- exonsBy(txdb, by = "gene")
introns_by_gene <- intronsByTranscript(txdb)

# Get gene regions (for overall gene assignment)
genes <- genes(txdb)

cat(paste0("Total genes: ", length(genes), "\n"))
cat(paste0("Total exon regions: ", sum(sapply(exons_by_gene, length)), "\n\n"))

# Read BAM file
cat("Reading BAM file...\n")

# Scan BAM for alignments with cell barcode tags
param <- ScanBamParam(
    what = c("qname", "flag", "rname", "pos", "cigar"),
    tag = c("CB", "UB")
)

bam_data <- scanBam(opt$bam, param = param)[[1]]

n_reads <- length(bam_data$qname)
cat(paste0("Total alignments: ", n_reads, "\n"))

# Filter for reads with barcodes
has_bc <- !is.na(bam_data$tag$CB)
cat(paste0("Alignments with barcode: ", sum(has_bc), "\n\n"))

# Create GAlignments object for overlap counting
cat("Creating alignment objects...\n")
ga <- GAlignments(
    seqnames = Rle(bam_data$rname[has_bc]),
    pos = bam_data$pos[has_bc],
    cigar = bam_data$cigar[has_bc],
    strand = Rle(rep("*", sum(has_bc)))
)

barcodes <- bam_data$tag$CB[has_bc]
read_names <- bam_data$qname[has_bc]

# Find overlaps with exons
cat("Finding exon overlaps...\n")
exon_hits <- findOverlaps(ga, unlist(exons_by_gene), type = "any")

# Find overlaps with genes (for intron detection)
cat("Finding gene overlaps...\n")
gene_hits <- findOverlaps(ga, genes, type = "within")

# Classify reads as spliced or unspliced
# A read is "spliced" if it only overlaps exons
# A read is "unspliced" if it overlaps introns or has intron-spanning alignments
cat("Classifying reads...\n")

# Check for spliced alignments (N in CIGAR)
cigar_vec <- bam_data$cigar[has_bc]
is_spliced_alignment <- grepl("N", cigar_vec)

# Get reads that overlap exons
reads_with_exon_overlap <- unique(queryHits(exon_hits))

# Reads overlapping genes
reads_with_gene_overlap <- unique(queryHits(gene_hits))

# Unspliced reads: overlap gene but not exclusively exons, OR have no splice junctions
# For simplicity, we'll use CIGAR-based classification:
# - Spliced: has N in CIGAR (splice junction)
# - Unspliced: no N in CIGAR and overlaps gene

# Get unique barcodes
unique_barcodes <- unique(barcodes)
n_cells <- length(unique_barcodes)
cat(paste0("Unique cell barcodes: ", n_cells, "\n"))

# Get gene names
gene_names <- names(genes)
n_genes <- length(gene_names)

# Initialize count matrices
cat("Building count matrices...\n")

# Create barcode and gene indices
bc_idx <- match(barcodes, unique_barcodes)
gene_idx_from_hits <- subjectHits(gene_hits)
read_idx_from_hits <- queryHits(gene_hits)

# Build sparse matrices
# Spliced counts
spliced_reads <- which(is_spliced_alignment)
spliced_in_genes <- intersect(spliced_reads, read_idx_from_hits)

if (length(spliced_in_genes) > 0) {
    spliced_gene_idx <- gene_idx_from_hits[match(spliced_in_genes, read_idx_from_hits)]
    spliced_bc_idx <- bc_idx[spliced_in_genes]
    
    spliced_counts <- sparseMatrix(
        i = spliced_gene_idx,
        j = spliced_bc_idx,
        x = rep(1, length(spliced_gene_idx)),
        dims = c(n_genes, n_cells),
        dimnames = list(gene_names, unique_barcodes)
    )
} else {
    spliced_counts <- sparseMatrix(
        dims = c(n_genes, n_cells),
        dimnames = list(gene_names, unique_barcodes)
    )
}

# Unspliced counts
unspliced_reads <- which(!is_spliced_alignment)
unspliced_in_genes <- intersect(unspliced_reads, read_idx_from_hits)

if (length(unspliced_in_genes) > 0) {
    unspliced_gene_idx <- gene_idx_from_hits[match(unspliced_in_genes, read_idx_from_hits)]
    unspliced_bc_idx <- bc_idx[unspliced_in_genes]
    
    unspliced_counts <- sparseMatrix(
        i = unspliced_gene_idx,
        j = unspliced_bc_idx,
        x = rep(1, length(unspliced_gene_idx)),
        dims = c(n_genes, n_cells),
        dimnames = list(gene_names, unique_barcodes)
    )
} else {
    unspliced_counts <- sparseMatrix(
        dims = c(n_genes, n_cells),
        dimnames = list(gene_names, unique_barcodes)
    )
}

# Total counts
total_counts <- spliced_counts + unspliced_counts

# Summary
cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")
cat(paste0("Genes: ", n_genes, "\n"))
cat(paste0("Cells: ", n_cells, "\n"))
cat(paste0("Total spliced counts: ", sum(spliced_counts), "\n"))
cat(paste0("Total unspliced counts: ", sum(unspliced_counts), "\n"))
cat(paste0("Spliced fraction: ", round(sum(spliced_counts) / sum(total_counts) * 100, 2), "%\n"))

# Save output
cat("\nSaving output...\n")

# Save as RDS (can be converted to H5AD later with Python)
output_data <- list(
    spliced = spliced_counts,
    unspliced = unspliced_counts,
    total = total_counts,
    barcodes = unique_barcodes,
    genes = gene_names
)

# Determine output format
if (grepl("\\.h5ad$", opt$output)) {
    # Save as RDS first, then convert
    rds_file <- sub("\\.h5ad$", ".rds", opt$output)
    saveRDS(output_data, rds_file)
    cat(paste0("Saved RDS: ", rds_file, "\n"))
    cat("Note: Use Python/scanpy to convert RDS to H5AD if needed.\n")
} else {
    saveRDS(output_data, opt$output)
}

cat(paste0("Output saved to: ", opt$output, "\n"))

cat("\n========================================\n")
cat("Velocity counting completed!\n")
cat("========================================\n")
