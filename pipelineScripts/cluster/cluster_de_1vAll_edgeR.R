library(openxlsx)
library(edgeR)
library(loomR)
library(parallel)

loom_file <- snakemake@input[["loom"]]
cluster_file <- snakemake@input[["clusters"]]

out_file <- snakemake@output[["out"]]
out_xlsx <- snakemake@output[["out_xlsx"]]

ds <- connect(filename = loom_file, mode = "r")

counts <- t(ds$matrix[, ])
barcodes <- ds$col.attrs$Barcode[]
ensIDs <- ds$row.attrs$EnsID[]

ds$close_all()

rownames(counts) <- ensIDs
colnames(counts) <- barcodes


clusters <- read.table(
    cluster_file, sep = "\t", header = TRUE, row.names = 1)

counts <- counts[, rownames(clusters)]

# Filter the genes
valid_genes <- rowSums(counts) > 0
counts <- counts[valid_genes, ]
counts <- data.matrix(counts)

de_cluster <- function(cluster_id){

    out_group <- setdiff(
        unique(clusters$Cluster),
        cluster_id
    )

    low_clusters <- out_group
    high_clusters <- cluster_id

    # Remove cells not in either comparison
    to_keep <- clusters$Cluster %in% union(low_clusters, high_clusters)
    counts_sub <- counts[, to_keep]
    clusters_sub <- clusters[to_keep, , drop = FALSE]

    # relabel groups as 'low' and 'high'
    group_labels <- ifelse(clusters_sub$Cluster %in% low_clusters,
                           "low", "high")
    group_labels <- as.factor(group_labels)
    group_labels <- relevel(group_labels, ref = "low")

    clusters_sub <- data.frame(group = group_labels,
                               row.names = rownames(clusters_sub))

    y <- DGEList(counts = counts_sub, samples = clusters_sub)
    y <- calcNormFactors(y, method = "none") # just use library size

    # Estimate dispersions
    design <- model.matrix(~clusters_sub$group)
    y <- estimateDisp(y, design) # 821 seconds for 7k cells

    fit <- glmFit(y, design)

    lrt <- glmLRT(fit) # Results stored in lrt$table

    results <- as.data.frame(lrt$table)
    results$FDR <- p.adjust(results$PValue, method = "BH")
    results <- results[order(results$FDR, -1 * abs(results$logFC)), ]

    results$cluster <- cluster_id

    results <- cbind(gene = rownames(results), results)
    rownames(results) <- NULL

    return(results)
}

all_cluster_ids <- unique(clusters$Cluster)

all_cluster_de <- mclapply(all_cluster_ids, de_cluster,
                           mc.cores = min(length(all_cluster_ids), 5)
                          )

all_concat <- do.call(rbind, all_cluster_de)

dir.create(dirname(out_file), showWarnings = FALSE)
write.table(all_concat,
            file = gzfile(out_file),
            row.names = FALSE,
            sep = "\t")

# Write results to excel also
wb <- createWorkbook()

for (cluster in unique(all_concat$cluster)){
    addWorksheet(wb, cluster)
    sub_data <- all_concat[all_concat$cluster == cluster, ]
    sub_data <- sub_data[order(sub_data$FDR, -1 * abs(sub_data$logFC)), ]
    writeData(wb, sheet = cluster, x = sub_data)
}

saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)
