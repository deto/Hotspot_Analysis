library(dynamicTreeCut)
X11.options(type = "Xlib")

results_file = "hotspot/hotspot_pairs_z.txt.gz"
results_hs_file = "hotspot/hotspot_hvg.txt"


dd <- read.table(gzfile(results_file), header=TRUE, sep="\t", row.names=1)
dd <- as.matrix(dd)

dd_dist <- dd * -1
offset <- min(dd_dist) * -1
dd_dist <- dd_dist + offset
diag(dd_dist) <- 0


dd_dist <- as.dist(dd_dist)

dendro <- hclust(dd_dist, method = "average")

MIN_CLUSTER_Z <- 3
maxHeight <- MIN_CLUSTER_Z * -1 + offset

clusters <- cutreeDynamicTree(
    dendro, maxTreeHeight = maxHeight,
    deepSplit = FALSE,
    minModuleSize = 5
    )

clusters <- as.factor(clusters)
unique(clusters)

library(NMF)

dd_ordered <- dd[dendro$order, dendro$order]
clusters_ordered <- clusters[dendro$order]

aheatmap(dd_ordered, Rowv = NA, Colv = NA,
    annRow = clusters_ordered,
    breaks = c(-1e99, seq(-25, 25, length.out = 99), 1e99)
)
