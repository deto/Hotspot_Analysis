library(WGCNA)
library(loomR)

# %% Load data

genes_file <- snakemake@input[["genes"]]

loom_file <- snakemake@input[["loom"]]

modules_out <- snakemake@output[["cluster_output"]]
scores_out <- snakemake@output[["scores"]]

minModuleSize <- snakemake@params[["min_module_size"]]

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
gene_ids <- ds$row.attrs$EnsID[]

ds$close_all()

# %% Subset to the thresholded genes

valid_genes <- readLines(genes_file)
# valid_genes <- readLines("genes/threshold.txt")
# valid_genes <- readLines("genes/hvg.txt")


is_valid_gene <- gene_ids %in% valid_genes

expression <- expression[is_valid_gene, , drop = FALSE]
gene_symbols <- gene_symbols[is_valid_gene]
gene_ids <- gene_ids[is_valid_gene]

rownames(expression) <- gene_ids
colnames(expression) <- rownames(meta)

datExpr <- log1p(t(expression))

# %% WGCNA Automatic tutorial
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-NetworkConstruction-auto.R


# %% =====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

if (snakemake@params[["power"]] == "auto"){

    # Choose a set of soft-thresholding powers
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    # Call the network topology analysis function
    sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
    # Plot the results:
    # sizeGrWindow(9, 5)
    # par(mfrow = c(1,2));
    # cex1 = 0.9;
    # # Scale-free topology fit index as a function of the soft-thresholding power
    # plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    #      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    #      main = paste("Scale independence"));
    # text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    #      labels=powers,cex=cex1,col="red");
    # # this line corresponds to using an R^2 cut-off of h
    # abline(h=0.90,col="red")
    # # Mean connectivity as a function of the soft-thresholding power
    # plot(sft$fitIndices[,1], sft$fitIndices[,5],
    #      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    #      main = paste("Mean connectivity"))
    # text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

    signedR2 <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
    power <- min(sft$fitIndices[signedR2 > .9, 1])
} else {
    power <- snakemake@params[["power"]]
}

message("Using power=", power)


# %% =====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

# Selected power = 3 based on the previous section's plots
net = blockwiseModules(datExpr, power = power,
                       TOMType = "unsigned", minModuleSize = minModuleSize,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3)


# %%=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# open a graphics window
# sizeGrWindow(12, 9)
# Convert labels to colors for plotting
# mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
# plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)


# %%=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
# save(MEs, moduleLabels, moduleColors, geneTree,
#      file = "FemaleLiver-02-networkConstruction-auto.RData")

modules <- data.frame(
    EnsID = names(moduleLabels),
    Cluster = moduleLabels
)
modules[modules$Cluster == 0, "Cluster"] <- -1

write.table(modules, modules_out, sep = "\t", row.names = FALSE)
write.table(MEs, scores_out, sep = "\t", col.names = NA)
