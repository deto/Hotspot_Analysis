library(fdrtool)
library(matrixStats)

# Load results

in_components_file <- snakemake@input[["components"]]
qvalcutoff <- snakemake@params[["qvalcutoff"]]
out_file <- snakemake@output[["components_fdr"]]
modules_file <- snakemake@output[["modules"]]

components <- read.table(
    in_components_file, sep = "\t", header = TRUE, row.names = 1
)

results <- lapply(seq(ncol(components)), function(i) {
    col_vals <- components[, i]
    col_result <- fdrtool(
        col_vals, plot = FALSE, cutoff.method = "fndr", verbose = FALSE
    )
    return(col_result$qval)
})


results <- do.call(cbind, results)

rownames(results) <- rownames(components)
colnames(results) <- colnames(components)

write.table(results, out_file, sep = "\t", col.names = NA)


# Use qvalcutoff to assign crisp modules
min_results <- rowMins(results)
results[results > qvalcutoff] <- 1e99
results[results > min_results] <- 1e99
results[results == min_results] <- 0

# Initialize to -1
modules <- numeric(nrow(results)) - 1

for (i in seq(ncol(results))) {
    modules[results[, i] == 0] <- i
}

modules <- data.frame(
    Cluster = modules, row.names = rownames(results)
)

write.table(modules, modules_file, sep = "\t", col.names = NA)
