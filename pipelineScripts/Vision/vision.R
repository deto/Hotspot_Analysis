devtools::load_all("/data/yosef2/users/david.detomaso/repos/VISION")
library(loomR)

loom_file <- snakemake@input[["loom"]]
out_file <- snakemake@output[["out"]]

out_dir <- dirname(out_file)
dir.create(out_dir, showWarnings = FALSE)

ds <- connect(filename = loom_file, mode = "r")

expression <- t(ds$layers$scaled[, ])

meta <- lapply(
    setNames(names(ds$col.attrs), names(ds$col.attrs)),
    function(col) {
        if (length(ds$col.attrs[[col]]$dims) > 1){
            NULL
        } else {
            ds$col.attrs[[col]][]
        }
    })
meta <- meta[!sapply(meta, is.null)]
meta <- data.frame(meta)
rownames(meta) <- meta$Barcode
meta$Barcode <- NULL
meta$LogNumUmi <- log10(meta$NumUmi)

gene_symbols <- ds$row.attrs$Symbol[]

ds$close_all()

expression <- convertGeneIds(expression, gene_symbols)
colnames(expression) <- rownames(meta)

sigs <- list()

if ("unnorm" %in% names(snakemake@input)){
    unnorm_file <- snakemake@input[["unnorm"]]
    if (endsWith(unnorm_file, ".feather")){
        unnorm <- read_feather_matrix(unnorm_file)
    } else {
        fr_input <- unnorm_file
        if (endsWith(fr_input, ".gz")){
            fr_input <- paste0("zcat ", fr_input)
        }
        unnorm <- fread(input = fr_input, sep = "\t",
                            header = TRUE, data.table = FALSE)

        rownames(unnorm) <- unnorm[, 1]
        unnorm[, 1] <- NULL
        unnorm <- as.matrix(unnorm)
    }
} else {
    unnorm <- NULL
}

if ("proteins" %in% names(snakemake@input)) {
    protein_file <- snakemake@input[["proteins"]]
    proteinData <- read.table(gzfile(protein_file), header = TRUE,
        row.names = 1)
    proteinData <- proteinData[colnames(expression), , drop = FALSE]
    proteinData <- log2(proteinData + 1)
} else {
    proteinData <- NULL
}

# Load any file that starts with 'meta_' as more metadata
for (input_name in names(snakemake@input)) {
    if (startsWith(input_name, "meta_")) {
        meta_file <- snakemake@input[[input_name]]
        newMeta <- read.delim(meta_file, check.names = FALSE, row.names = 1)
        meta <- merge(meta, newMeta, by = "row.names")
        row.names(meta) <- meta$Row.names
        meta$Row.names <- NULL
    }
}

# Load latent space coordinates
if ("latent" %in% names(snakemake@input)){
    latent_file <- snakemake@input[["latent"]]
    latentSpace <- read.delim(latent_file, check.names = FALSE, row.names = 1)
} else {
    latentSpace <- NULL
}

# Name the results
wd <- getwd()
name <- basename(wd)

projection_methods <- c("tSNE30", "UMAP")
if ("tsne" %in% names(snakemake@input) || "umap" %in% names(snakemake@input) ||
    any(startsWith(names(snakemake@input), "proj_"))) {
    projection_methods <- character()
}

fp <- Vision(expression, signatures = sigs, meta = meta,
             projection_genes = rownames(expression),
             proteinData = proteinData,
             name = name,
             unnormalizedData = unnorm,
             projection_methods = projection_methods,
             latentSpace = latentSpace, pool = FALSE)

# Load projection coordinates
if ("tsne" %in% names(snakemake@input)) {
    proj_file <- snakemake@input[["tsne"]]
    tsne <- read.delim(proj_file, check.names = FALSE, row.names = 1)
    fp <- addProjection(fp, "TSNE", tsne)
}

if ("umap" %in% names(snakemake@input)) {
    proj_file <- snakemake@input[["umap"]]
    umap <- read.delim(proj_file, check.names = FALSE, row.names = 1)
    fp <- addProjection(fp, "UMAP", umap)
}

# Load any file that starts with 'proj_' as a projection
for (input_name in names(snakemake@input)) {
    if (startsWith(input_name, "proj_")) {
        proj_name <- gsub("proj_", "", input_name)
        proj_file <- snakemake@input[[input_name]]
        coords <- read.delim(proj_file, check.names = FALSE, row.names = 1)
        fp <- addProjection(fp, proj_name, coords)
    }
}

options(mc.cores = 10)
fp <- analyze(fp)

saveRDS(fp, out_file)
