library(Seurat)
library(feather)


counts_file <- snakemake@input[["expression"]]
# counts_file <- "../../data/CiteSeq_PBMC/expression_counts.feather"

genes_out_file <- snakemake@output[["genes"]]

counts <- feather::read_feather(counts_file)

counts <- as.data.frame(counts)
rownames(counts) <- counts[["index"]]
counts$index <- NULL
counts <- data.matrix(counts)


if ("doXcutoff" %in% names(snakemake@params) &&
    snakemake@params[["doXcutoff"]]) {

    x.cutoff <- 3
} else {
    x.cutoff <- 100
}

obj <- CreateSeuratObject(raw.data = counts, min.cells = 10, min.features = 0)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

obj <- FindVariableGenes(object = obj, mean.function = ExpMean, dispersion.function = LogVMR,
    x.low.cutoff = 0.0125, x.high.cutoff = x.cutoff, y.cutoff = 0.5, do.plot = FALSE)

genes <- obj@var.genes


writeLines(genes, file(genes_out_file))
