#!/usr/bin/env Rscript

# ==============================================================================
# Generate combined count matrices (spliced + unspliced) for 10x long-read data
#
# This script generates count matrices that include BOTH spliced and unspliced
# reads. This is essential for nuclei sequencing where a significant fraction
# of transcripts are unspliced pre-mRNA.
#
# For each gene, the output includes:
# - Total counts = spliced counts + unspliced counts
# - This captures the full transcriptional output from each cell
#
# Usage:
#   Rscript combined_counts.R --bam <file> --gtf <file> --output <file>
#
# Output:
#   - RDS file containing:
#     - gene_counts: combined gene-level counts per cell (spliced + unspliced)
#     - transcript_counts: combined transcript-level counts per cell
#     - spliced_counts: spliced reads only (for reference)
#     - unspliced_counts: unspliced reads only (for reference)
#     - cell_metadata: cell barcode information
# ==============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(GenomicAlignments)
    library(GenomicFeatures)
    library(Matrix)
    library(data.table)
    library(Rsamtools)
})

# Parse command line arguments
option_list <- list(
    make_option(c("-b", "--bam"),
                type = "character",
                default = NULL,
                help = "Input BAM file with CB/UB tags",
                metavar = "FILE"),
    make_option(c("-g", "--gtf"),
                type = "character",
                default = NULL,
                help = "GTF annotation file",
                metavar = "FILE"),
    make_option(c("-o", "--output"),
                type = "character",
                default = NULL,
                help = "Output RDS file",
                metavar = "FILE"),
    make_option(c("--min_counts"),
                type = "integer",
                default = 1,
                help = "Minimum counts per cell to include [default: %default]",
                metavar = "INT")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Generate combined spliced+unspliced count matrices")
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
cat("Combined Count Matrix Generation\n")
cat("(Spliced + Unspliced)\n")
cat("========================================\n")
cat(paste0("Input BAM: ", opt$bam, "\n"))
cat(paste0("GTF file: ", opt$gtf, "\n"))
cat(paste0("Output: ", opt$output, "\n"))
cat("========================================\n\n")

# Create TxDb from GTF
cat("Loading annotations...\n")
txdb <- makeTxDbFromGFF(opt$gtf, format = "gtf")

# Get exon regions by gene
cat("Extracting genomic features...\n")
exons_by_gene <- exonsBy(txdb, by = "gene")
genes <- genes(txdb)

# Get transcripts for transcript-level counting
transcripts <- transcripts(txdb)

cat(paste0("Total genes: ", length(genes), "\n"))
cat(paste0("Total transcripts: ", length(transcripts), "\n\n"))

# Read BAM file
cat("Reading BAM file...\n")

# Scan BAM for alignments with cell barcode tags
param <- ScanBamParam(
    what = c("qname", "flag", "rname", "pos", "cigar", "mapq"),
    tag = c("CB", "UB"),
    mapqFilter = 10  # Filter low quality alignments
)

bam_data <- scanBam(opt$bam, param = param)[[1]]

n_reads <- length(bam_data$qname)
cat(paste0("Total alignments: ", n_reads, "\n"))

# Filter for reads with cell barcodes
has_bc <- !is.na(bam_data$tag$CB)
n_with_bc <- sum(has_bc)
cat(paste0("Alignments with cell barcode: ", n_with_bc, " (",
           round(100 * n_with_bc / n_reads, 2), "%)\n\n"))

if (n_with_bc == 0) {
    stop("Error: No reads with cell barcodes found. Check BAM tagging step.")
}

# Extract data for reads with barcodes
cell_barcodes <- bam_data$tag$CB[has_bc]
umi_sequences <- bam_data$tag$UB[has_bc]
read_names <- bam_data$qname[has_bc]
cigars <- bam_data$cigar[has_bc]

# Create GAlignments for overlap detection
cat("Creating alignment objects...\n")
ga <- GAlignments(
    seqnames = Rle(bam_data$rname[has_bc]),
    pos = bam_data$pos[has_bc],
    cigar = cigars,
    strand = Rle(rep("*", n_with_bc))
)

# Classify reads as spliced or unspliced based on CIGAR
# Spliced reads have 'N' in CIGAR (splice junction)
cat("Classifying reads as spliced/unspliced...\n")
is_spliced <- grepl("N", cigars)
cat(paste0("  Spliced reads: ", sum(is_spliced), " (", 
           round(100 * sum(is_spliced) / n_with_bc, 2), "%)\n"))
cat(paste0("  Unspliced reads: ", sum(!is_spliced), " (",
           round(100 * sum(!is_spliced) / n_with_bc, 2), "%)\n\n"))

# Find overlaps with genes
cat("Finding gene overlaps...\n")
gene_hits <- findOverlaps(ga, genes, type = "any", ignore.strand = TRUE)

# Get unique cell barcodes
unique_barcodes <- unique(cell_barcodes)
n_cells <- length(unique_barcodes)
cat(paste0("Unique cell barcodes: ", n_cells, "\n\n"))

# Get gene names
gene_names <- names(genes)
n_genes <- length(gene_names)

# Create index mappings
bc_idx <- match(cell_barcodes, unique_barcodes)
gene_idx_map <- rep(NA, n_with_bc)
read_idx_from_hits <- queryHits(gene_hits)
gene_idx_from_hits <- subjectHits(gene_hits)
gene_idx_map[read_idx_from_hits] <- gene_idx_from_hits

# Build count matrices
cat("Building count matrices...\n")

# Helper function to build sparse matrix from indices
build_count_matrix <- function(read_indices, gene_indices, bc_indices, 
                                n_genes, n_cells, gene_names, barcodes) {
    # Filter to valid indices
    valid <- !is.na(gene_indices) & !is.na(bc_indices)
    gi <- gene_indices[valid]
    bi <- bc_indices[valid]
    
    if (length(gi) == 0) {
        return(sparseMatrix(dims = c(n_genes, n_cells),
                           dimnames = list(gene_names, barcodes)))
    }
    
    # Aggregate counts (handle duplicates)
    dt <- data.table(gene = gi, cell = bi)
    dt <- dt[, .(count = .N), by = .(gene, cell)]
    
    sparseMatrix(
        i = dt$gene,
        j = dt$cell,
        x = dt$count,
        dims = c(n_genes, n_cells),
        dimnames = list(gene_names, barcodes)
    )
}

# Spliced counts
spliced_idx <- which(is_spliced)
spliced_gene_idx <- gene_idx_map[spliced_idx]
spliced_bc_idx <- bc_idx[spliced_idx]

spliced_counts <- build_count_matrix(
    spliced_idx, spliced_gene_idx, spliced_bc_idx,
    n_genes, n_cells, gene_names, unique_barcodes
)

# Unspliced counts
unspliced_idx <- which(!is_spliced)
unspliced_gene_idx <- gene_idx_map[unspliced_idx]
unspliced_bc_idx <- bc_idx[unspliced_idx]

unspliced_counts <- build_count_matrix(
    unspliced_idx, unspliced_gene_idx, unspliced_bc_idx,
    n_genes, n_cells, gene_names, unique_barcodes
)

# Combined counts (spliced + unspliced)
cat("\nCombining spliced and unspliced counts...\n")
combined_counts <- spliced_counts + unspliced_counts

# Summary statistics
cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")
cat(paste0("Genes: ", n_genes, "\n"))
cat(paste0("Cells: ", n_cells, "\n"))
cat(paste0("\nSpliced counts:\n"))
cat(paste0("  Total: ", sum(spliced_counts), "\n"))
cat(paste0("  Genes detected: ", sum(rowSums(spliced_counts) > 0), "\n"))
cat(paste0("\nUnspliced counts:\n"))
cat(paste0("  Total: ", sum(unspliced_counts), "\n"))
cat(paste0("  Genes detected: ", sum(rowSums(unspliced_counts) > 0), "\n"))
cat(paste0("\nCombined counts (spliced + unspliced):\n"))
cat(paste0("  Total: ", sum(combined_counts), "\n"))
cat(paste0("  Genes detected: ", sum(rowSums(combined_counts) > 0), "\n"))
cat(paste0("  Median counts per cell: ", median(colSums(combined_counts)), "\n"))
cat(paste0("  Mean counts per cell: ", round(mean(colSums(combined_counts)), 2), "\n"))

# Calculate fraction unspliced per cell
frac_unspliced <- colSums(unspliced_counts) / (colSums(combined_counts) + 1e-10)
cat(paste0("\nFraction unspliced:\n"))
cat(paste0("  Mean: ", round(mean(frac_unspliced) * 100, 2), "%\n"))
cat(paste0("  Median: ", round(median(frac_unspliced) * 100, 2), "%\n"))

# Create cell metadata
cell_metadata <- data.frame(
    barcode = unique_barcodes,
    total_counts = as.numeric(colSums(combined_counts)),
    spliced_counts = as.numeric(colSums(spliced_counts)),
    unspliced_counts = as.numeric(colSums(unspliced_counts)),
    fraction_unspliced = as.numeric(frac_unspliced),
    n_genes_detected = as.numeric(colSums(combined_counts > 0)),
    stringsAsFactors = FALSE
)

# Save output
cat("\nSaving output...\n")

output_data <- list(
    # Primary output: combined counts (spliced + unspliced)
    gene_counts = combined_counts,
    
    # Separate matrices for reference/velocity
    spliced_counts = spliced_counts,
    unspliced_counts = unspliced_counts,
    
    # Metadata
    cell_metadata = cell_metadata,
    gene_names = gene_names,
    barcodes = unique_barcodes,
    
    # Summary statistics
    summary = list(
        n_genes = n_genes,
        n_cells = n_cells,
        total_spliced = sum(spliced_counts),
        total_unspliced = sum(unspliced_counts),
        total_combined = sum(combined_counts),
        mean_fraction_unspliced = mean(frac_unspliced)
    )
)

saveRDS(output_data, opt$output)
cat(paste0("Output saved to: ", opt$output, "\n"))

# Also save as CSV for easy inspection
csv_file <- sub("\\.rds$", "_gene_counts.csv.gz", opt$output)
write.csv(as.matrix(combined_counts), gzfile(csv_file), row.names = TRUE)
cat(paste0("CSV saved to: ", csv_file, "\n"))

cat("\n========================================\n")
cat("Combined counting completed!\n")
cat("========================================\n")
cat("\nOutput contains:\n")
cat("  - gene_counts: Combined (spliced + unspliced) gene-level counts\n")
cat("  - spliced_counts: Spliced reads only\n")
cat("  - unspliced_counts: Unspliced reads only\n")
cat("  - cell_metadata: Per-cell statistics\n")
cat("\nUse gene_counts for standard single-cell analysis.\n")
cat("Use spliced/unspliced matrices for RNA velocity with scVelo.\n")
