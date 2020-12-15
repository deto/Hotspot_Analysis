devtools::load_all("/data/yosef2/users/david.detomaso/repos/SymSim/")
library(Rtsne)
library(loomR)
library(ggplot2)
library(pbmcapply)

true_counts_file <- snakemake@input[["true_counts"]]
cell_meta_file <- snakemake@input[["cell_meta"]]


depth_mean <- as.numeric(snakemake@params[["seq_depth"]])
depth_sd <- depth_mean / 10

obs_counts_file <- snakemake@output[["obs_counts"]]


true_counts <- read.table(gzfile(true_counts_file), sep = "\t",
    row.names = 1, header = TRUE)
true_counts <- as.matrix(true_counts)


cell_meta <- read.table(cell_meta_file, sep = "\t", header = TRUE)


# What about the technical counts then

ngenes_total <- nrow(true_counts)

path <- system.file("data/gene_len_pool.RData", package = "SymSim")
load(path)
gene_len <- sample(gene_len_pool, ngenes_total, replace = FALSE)

options(mc.cores = 20)
observed_counts_res <- True2ObservedCounts(
    true_counts = true_counts,
    meta_cell = cell_meta,
    protocol = "UMI",
    alpha_mean = 0.05, alpha_sd = 0.02, lenslope = 0.01,
    nbins = 20, gene_len = gene_len,
    amp_bias_limit = c(-0.3, 0.3), rate_2PCR = 0.7,
    nPCR = 16, depth_mean = depth_mean, depth_sd = depth_sd
)

observed_counts <- observed_counts_res[[1]]


# %% create output loom file

# Loom files broken in R, make it in python instead...
observed_counts_df <- data.frame(observed_counts)
rownames(observed_counts_df) <- rownames(true_counts)
colnames(observed_counts_df) <- colnames(true_counts)

write.table(
    observed_counts_df, gzfile(obs_counts_file),
    sep = "\t", col.names = NA)
