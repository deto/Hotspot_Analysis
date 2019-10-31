library(M3Drop)
library(loomR)

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

count_mat <- NBumiConvertData(counts, is.counts = TRUE)

DANB_fit <- NBumiFitModel(count_mat)

NBDropFS <- NBumiFeatureSelectionCombinedDrop(
    DANB_fit, method = "fdr",
    ntop = nrow(count_mat),
)

genes <- as.character(
    NBDropFS[NBDropFS$q.value < 0.05, "Gene"]
)


if ("highXMeanCutoff" %in% names(snakemake@params)) {
    high_cut <- snakemake@params[["highXMeanCutoff"]]
    scaled_counts <- t(t(counts) / colSums(counts)) * 10000
    geneMeans <- rowMeans(scaled_counts)
    valid_genes <- names(geneMeans)[geneMeans < high_cut]
    genes <- intersect(genes, valid_genes)
}

if ("geneInfo" %in% names(snakemake@output)) {
    geneInfoFile <- snakemake@output[["geneInfo"]]
    write.table(NBDropFS, geneInfoFile, sep = "\t")
}

writeLines(genes, file(genes_out_file))
