library(Seurat)
library(loomR)


loom_file <- snakemake@input[["loom"]]
genes_out_file <- snakemake@output[["genes"]]

ds <- connect(filename = loom_file, mode = "r")

counts <- t(ds$matrix[, ])
barcodes <- ds$col.attrs$Barcode[]
gene_symbols <- ds$row.attrs$EnsID[]

ds$close_all()

rownames(counts) <- gene_symbols
colnames(counts) <- barcodes

if ("doXcutoff" %in% names(snakemake@params) &&
    snakemake@params[["doXcutoff"]]) {
    x.cutoff <- 3
} else {
    x.cutoff <- 100
}

if ("lowXcutoff" %in% names(snakemake@params)) {
    lowXcutoff <- snakemake@params[["lowXcutoff"]]
} else {
    lowXcutoff <- .0125
}

obj <- CreateSeuratObject(raw.data = counts, min.cells = 10, min.features = 0)

obj <- NormalizeData(
    obj, normalization.method = "LogNormalize",
    scale.factor = 10000)

obj <- FindVariableGenes(object = obj, mean.function = ExpMean,
    dispersion.function = LogVMR, x.low.cutoff = lowXcutoff,
    x.high.cutoff = x.cutoff, y.cutoff = 0.5, do.plot = FALSE)

genes <- obj@var.genes

if ("highXMeanCutoff" %in% names(snakemake@params)) {
    high_cut <- snakemake@params[["highXMeanCutoff"]]
    geneMeans <- rowMeans(exp(obj@data) - 1)
    valid_genes <- names(geneMeans)[geneMeans < high_cut]
    genes <- intersect(genes, valid_genes)
}

if ("geneInfo" %in% names(snakemake@output)) {
    geneInfoFile <- snakemake@output[["geneInfo"]]
    write.table(obj@hvg.info, geneInfoFile, sep = "\t")
}

writeLines(genes, file(genes_out_file))
