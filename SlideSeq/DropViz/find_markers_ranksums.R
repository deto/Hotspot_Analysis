devtools::load_all("DropSeq.util")
devtools::load_all("/data/yosef2/users/david.detomaso/repos/VISION")


dge <- loadSparseDge("F_GRCm38.81.P60Cerebellum_ALT.raw.dge.txt.gz")

ca <- readRDS("F_GRCm38.81.P60Cerebellum_ALT.cluster.assign.RDS")

dge <- dge[, names(ca)]

valid_genes <- rowSums(dge) > 0
dge <- dge[valid_genes, ]

ca <- as.data.frame(ca)
dge <- as.matrix(dge)


options(mc.cores = 16)

rr <- lapply(setNames(levels(ca$ca), levels(ca$ca)),
    function(level){

    ca$Group <- "Denom"
    ca$Group[ca$ca == level] <- "Num"
    ca$Group <- relevel(factor(ca$Group), ref = "Denom")
    Group <- ca$Group


    cluster_num <- which(ca$Group == "Num")
    cluster_denom <- which(ca$Group == "Denom")

    res <- matrix_wilcox_cpp(dge, cluster_num, cluster_denom)

    res$FDR <- p.adjust(res$pval, method = "BH")

    num_mean <- (rowSums(dge[, cluster_num]) + 1) / length(cluster_num)
    denom_mean <- (rowSums(dge[, cluster_denom]) + 1) / length(cluster_denom)
    res$logFC <- log2(num_mean / denom_mean)

    res$level <- level
    return(res)
})


all_res_cat <- do.call(rbind, rr)
write.table(all_res_cat, "ranksums_markers_1vAll.txt", sep = "\t")


topGenes <- function(mat){
    msub <- mat[mat$logFC > 0, ]
    msub <- msub[order(msub$AUC), ]
    return(msub)
}
