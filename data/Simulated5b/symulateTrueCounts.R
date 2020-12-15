devtools::load_all("/data/yosef2/users/david.detomaso/repos/SymSim/")
library(Rtsne)
library(loomR)
library(ggplot2)
library(pbmcapply)
library(phytools)

gene_effects_file <- snakemake@output[["gene_effects"]]
true_counts_file <- snakemake@output[["true_counts"]]
cell_meta_file <- snakemake@output[["cell_meta"]]

# %% Follow the vignette

phyla <- read.tree(
    text = "(((A:1,B:1):1,(C:0.5,D:0.5):1.5):1,E:3);"
)

true_counts_res <- SimulateTrueCounts(
    ncells_total = 3000,
    min_popsize = 200,
    ngenes = 5000,
    evf_center = 1,
    evf_type = "discrete",
    nevf = 10,
    n_de_evf = 5,
    Sigma = 0.5,
    vary = "s",
    gene_effect_prob = 0.3,
    geffect_mean = 0,
    gene_effects_sd = 1,
    param_realdata = "zeisel.imputed",
    phyla = phyla,
    randseed = 8675309
)

true_counts <- true_counts_res[[1]]
gene_effects <- true_counts_res[[2]]
cell_meta <- true_counts_res[[3]]
params <- true_counts_res[[4]]



# %% create output files

gene_ids <- paste0("Gene", seq(nrow(true_counts)))
cell_ids <- cell_meta$cellid

# True counts

true_counts_df <- data.frame(true_counts)
rownames(true_counts_df) <- gene_ids
colnames(true_counts_df) <- cell_ids

write.table(
    true_counts_df, gzfile(true_counts_file), sep = "\t", col.names = NA
)

# Component effects (for S)
gene_effects_df <- data.frame(gene_effects[[3]])
rownames(gene_effects_df) <- gene_ids
write.table(gene_effects_df, gene_effects_file, sep = "\t", col.names = NA)

# Save all the component effects this time!
gene_effects_df <- data.frame(gene_effects[[1]])
rownames(gene_effects_df) <- gene_ids
write.table(gene_effects_df, "true_counts/gene_effects_kon.txt",
    sep = "\t", col.names = NA
)

gene_effects_df <- data.frame(gene_effects[[2]])
rownames(gene_effects_df) <- gene_ids
write.table(gene_effects_df, "true_counts/gene_effects_koff.txt",
    sep = "\t", col.names = NA
)

gene_effects_df <- data.frame(gene_effects[[3]])
rownames(gene_effects_df) <- gene_ids
write.table(gene_effects_df, "true_counts/gene_effects_s.txt",
    sep = "\t", col.names = NA
)

# cell meta -> need this to generate downsampled reads
write.table(cell_meta, cell_meta_file, sep = "\t", row.names = FALSE)
