# Performs PCA-based feature selection
# Uses the implementation in M3Drop

library(M3Drop)
library(loomR)
library(matrixStats)

# loom_file <- "../../data/10x_PBMC_w_proteins/cd4/data.loom"

loom_file <- snakemake@input[["loom"]]
genes_out_file <- snakemake@output[["genes"]]


ds <- connect(filename = loom_file, mode = "r")

counts <- t(ds$matrix[, ])
barcodes <- ds$col.attrs$Barcode[]
gene_symbols <- ds$row.attrs$EnsID[]

ds$close_all()

rownames(counts) <- gene_symbols
colnames(counts) <- barcodes

scaled_counts <- t(t(counts) / colSums(counts)) * 10000

# Remove genes that are constant across cells
scaled_counts <- scaled_counts[
    rowVars(scaled_counts) > 0,
]

results <- irlbaPcaFS(scaled_counts, pcs = seq(5))

results <- results * -1  # Under their inversion

results <- data.frame(
    Score = results
)
results <- results[order(results$Score * -1), , drop = FALSE]


genes <- rownames(head(results, 500))

if ("highXMeanCutoff" %in% names(snakemake@params)) {
    high_cut <- snakemake@params[["highXMeanCutoff"]]
    geneMeans <- rowMeans(scaled_counts)
    valid_genes <- names(geneMeans)[geneMeans < high_cut]
    genes <- intersect(genes, valid_genes)
}

if ("geneInfo" %in% names(snakemake@output)) {
    geneInfoFile <- snakemake@output[["geneInfo"]]
    write.table(results, geneInfoFile, sep = "\t")
}

writeLines(genes, file(genes_out_file))
