devtools::load_all("DropSeq.util")
library(edgeR)


dge <- loadSparseDge("F_GRCm38.81.P60Cerebellum_ALT.raw.dge.txt.gz")

ca <- readRDS("F_GRCm38.81.P60Cerebellum_ALT.subcluster.assign.RDS")

dge <- dge[, names(ca)]

valid_genes <- rowSums(dge) > 0
dge <- dge[valid_genes, ]

ca <- as.data.frame(ca)

y <- DGEList(counts = as.matrix(dge), samples = ca)
y <- calcNormFactors(y, method = "none") # just use library size

options(mc.cores = 16)

# Ok, so it's clear that this doesn't work well.  Instead, we should do 1 vs. all comparisons and just take the top N genes per group as our set

all_res <- lapply(levels(ca$ca), function(level){

    print(paste0(level, Sys.time()))

    ca$Group <- "Denom"
    ca$Group[ca$ca == level] <- "Num"
    ca$Group <- relevel(factor(ca$Group), ref = "Denom")
    Group <- ca$Group

    design <- model.matrix(~Group)

    # 821 seconds for 7k cells
    # 1990 seconds for 12.5k cells, 21.4k genes
    y <- estimateDisp(y, design)

    fit <- glmFit(y, design) # 133 seconds for 12.5k cells, 21.4k genes

    lrt <- glmLRT(fit, coef = "GroupNum") # 75 seconds for 12.5k cells, 21.4k genes


    results <- as.data.frame(lrt$table)
    results$FDR <- p.adjust(results$PValue, method = "BH")
    results <- results[order(results$FDR), ]
    results$Cluster <- level
    GeneSymbols <- data.frame(GeneSymbol = row.names(results))
    results <- cbind(GeneSymbols, results)
    row.names(results) <- NULL
    return(results)
})

all_res_cat <- do.call(rbind, all_res)
write.table(all_res_cat, "edger_markers_1vAll_subcluster.txt", sep = "\t")
