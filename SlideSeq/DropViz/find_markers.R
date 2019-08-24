devtools::load_all("DropSeq.util")
library(edgeR)


dge <- loadSparseDge("F_GRCm38.81.P60Cerebellum_ALT.raw.dge.txt.gz")

ca <- readRDS("F_GRCm38.81.P60Cerebellum_ALT.cluster.assign.RDS")

dge <- dge[, names(ca)]

valid_genes <- rowSums(dge) > 0
dge <- dge[valid_genes, ]

ca <- as.data.frame(ca)

y <- DGEList(counts = as.matrix(dge), samples = ca)
y <- calcNormFactors(y, method = "none") # just use library size

options(mc.cores=16)

design <- model.matrix(~ca$ca)
y <- estimateDisp(y, design) # 821 seconds for 7k cells

fit <- glmQLFit(y, design)

coef <- c(
    "ca$ca2", "ca$ca3", "ca$ca4", "ca$ca5",
    "ca$ca6", "ca$ca7", "ca$ca8", "ca$ca9",
    "ca$ca10", "ca$ca11"
)

lrt <- glmQLFTest(fit, coef=coef)


results <- as.data.frame(lrt$table)
results$FDR <- p.adjust(results$PValue, method = "BH")
results <- results[order(results$FDR), ]

write.table(results, "edger_markers.txt", sep = "\t")



glmfit <- glmFit(y, design)

coef <- c(
    "ca$ca2", "ca$ca3", "ca$ca4", "ca$ca5",
    "ca$ca6", "ca$ca7", "ca$ca8", "ca$ca9",
    "ca$ca10", "ca$ca11"
)

lrt <- glmLRT(glmfit, coef=coef)


results_lrt <- as.data.frame(lrt$table)
results_lrt$FDR <- p.adjust(results_lrt$PValue, method = "BH")
results_lrt <- results_lrt[order(results_lrt$FDR), ]

write.table(results_lrt, "edger_markers_lrt.txt", sep = "\t")

# Ok, so it's clear that this doesn't work well.  Instead, we should do 1 vs. all comparisons and just take the top N genes per group as our set

all_res <- lapply(levels(ca$ca), function(level){

    ca$Group <- "Denom"
    ca$Group[ca$ca == level] <- "Num"
    ca$Group <- relevel(factor(ca$Group), ref = "Denom")
    Group <- ca$Group

    design <- model.matrix(~Group)
    y <- estimateDisp(y, design) # 821 seconds for 7k cells

    fit <- glmFit(y, design)

    lrt <- glmLRT(fit, coef = "GroupNum")


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
wite.table(all_res_cat, "edger_markers_1vAll.txt", sep = "\t")
