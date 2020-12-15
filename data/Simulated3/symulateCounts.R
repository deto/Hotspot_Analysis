devtools::load_all("/data/yosef2/users/david.detomaso/repos/SymSim/")
library(Rtsne)
library(loomR)
library(ggplot2)
library(pbmcapply)

cell_effects_file <- snakemake@output[["cell_effects"]]
gene_effects_file <- snakemake@output[["gene_effects"]]
gene_indices_file <- snakemake@output[["gene_indices"]]
true_counts_file <- snakemake@output[["true_counts"]]
obs_counts_file <- snakemake@output[["obs_counts"]]

# Inside this file I've broken the SimulateTrueCounts function into three pieces
# First piece is 'compute_evf_params': unchanged from main code
#    This results in evf_res and gene_effects
# Then I modify these data for my experiment
# Then trueCounts is the rest of SimulateTrueCounts

# library(SymSim)
# What are the variables of the simulation that we need to control?
# Number of Cells
# Number of Genes
# Number of Genes involved in confounding signal
# Overlap of genes in confounding signal and other EVF genes

# Compute evf_res and gene_effects
compute_evf_params <- function(ncells_total, ngenes, geffect_mean = 0,
                               gene_effects_sd = 1, gene_effect_prob = 0.3,
                               evf_center = 1, Sigma = 0.5, nevf = 10,
                               min_popsize = 200, randseed = 1, i_minpop = 1,
                               evf_type = "one.population", n_de_evf = 0, impulse = F,
                               vary = "all", phyla = NULL
                               ){

    set.seed(randseed)
    n_nd_evf = nevf - n_de_evf
    seed <- sample(c(1:1e5),size = 2)
    param_names <- c("kon", "koff", "s")

    if (evf_type == 'one.population'){
        evf_mean=rep(evf_center,nevf); evf_sd=rep(Sigma,nevf)
        evfs <- lapply(1:3, function(iparam){
            evfs_per_param <- lapply(c(1:ncells_total),function(celli){
                evf <- sapply(c(1:nevf),function(ievf){rnorm(1,evf_mean[ievf],evf_sd[ievf])})
                return(evf)})
            evfs_per_param <- do.call(rbind, evfs_per_param)
            colnames(evfs_per_param) <- sprintf("%s_evf%d", param_names[iparam], 1:nevf)
            return(evfs_per_param)})
        evf_res <- list(evfs=evfs, meta=data.frame(pop=rep(1, ncells_total)))
    } else if(evf_type=='discrete'){

        evf_res <- DiscreteEVF(phyla,ncells_total,min_popsize,i_minpop=i_minpop,Sigma,
                               n_nd_evf, n_de_evf, vary=vary, evf_center=evf_center, seed=seed[1])

    }else if(evf_type=='continuous'){

        evf_res <- ContinuousEVF(phyla,ncells_total,n_nd_evf=nevf-n_de_evf,n_de_evf=n_de_evf,
                                 evf_center=evf_center,vary=vary,impulse=impulse,
                                 Sigma,plotting=T,seed=seed[1])    
    }

    gene_effects <- GeneEffects(ngenes=ngenes, nevf=nevf, randseed=seed[2], prob=gene_effect_prob,
                                geffect_mean=geffect_mean, geffect_sd=gene_effects_sd)

    return(list(evf_res, gene_effects))
     
}

trueCounts <- function(evf_res, gene_effects) {

    bimod <- 0
    param_realdata <- "zeisel.imputed"
    ngenes <- nrow(gene_effects[[1]])
    ncells_total <- nrow(evf_res[[1]][[1]])
    joint <- F

    if (!is.null(param_realdata)){
        if (param_realdata == "zeisel.imputed"){
            path <- system.file("data/param_realdata.zeisel.imputed.robj",
                                package = "SymSim")
            load(path)
        }
        if (param_realdata == "zeisel.pop4"){
            path <- system.file("data/param_realdata.zeisel.pop4.robj",
                                package = "SymSim")
            load(path)
        }
        match_params[, 1] <- log(base = 10, match_params[, 1])
        match_params[, 2] <- log(base = 10, match_params[, 2])
        match_params_den <- lapply(c(1:3), function(i){
            density(match_params[, i], n = 2000)
        })
        if (joint == F){
            params <- Get_params(gene_effects, evf_res[[1]],
                                 match_params_den, bimod)
        } else {
            params <- Get_params3(gene_effects, evf_res[[1]],
                                  match_params, bimod)
        }
    } else {
        params <- Get_params2(gene_effects, evf_res[[1]], bimod, ranges)
    }

    counts <- pbmcapply::pbmclapply(c(1:ngenes), function(i){
        count <- sapply(c(1:ncells_total), function(j){
            y <- rbeta(1, params[[1]][i, j], params[[2]][i, j])
            x <- rpois(1, y * params[[3]][i, j])
            return(x)
        })
    })

    cell_meta <- cbind(cellid = paste("cell", seq(1, ncells_total), sep = "_"),
                       evf_res[[2]], evf_res[[1]])
    counts <- do.call(rbind, counts)
    result <- list(counts, gene_effects, cell_meta, params)
    return(result)
}

# %% Great - ok, now see if we can create multiple components

ncells_total <- 3000
ngenes_total <- 5000
nevf <- 1
min_popsize = 200
randseed = 1
i_minpop = 1
evf_type = "one.population"
n_de_evf = 0
impulse = F
vary = "s"
phyla = NULL
geffect_mean_0 <- 0
geffect_sd_0 <- 1  
evf_center <- 1
Sigma <- 0.5

rr <- compute_evf_params(
        ncells_total,
        ngenes = ngenes_total,
        geffect_mean = geffect_mean_0, gene_effects_sd = geffect_sd_0,
        evf_center = evf_center, Sigma = Sigma, nevf = nevf,
        min_popsize = min_popsize, randseed = randseed, i_minpop = i_minpop,
        evf_type = evf_type, n_de_evf = n_de_evf, impulse = impulse,
        vary = vary, phyla = phyla
)


evf_res <- rr[[1]]
evfs <- evf_res[[1]]
meta <- evf_res[[2]]
gene_effects <- rr[[2]]

# The first EVF should just be an offset term with random gene contributions
for (p in seq(3)){
    gene_effects[[p]][, 1] <- rnorm(nrow(gene_effects[[p]]))*geffect_sd_0 + geffect_mean_0
    evfs[[p]][, 1] <- rep(1, nrow(evfs[[p]]))
}

# Effect sizes for the new effects
new_effect_mean <- 0

strength <- as.numeric(snakemake@params[["strength"]])
new_effect_sd <- .1 * strength

n_evfs_to_add <- 5
n_genes_per_evf <- 100
n_cells_per_evf <- 500

gene_indices <- sample.int(
    ngenes_total,
    size = n_evfs_to_add * n_genes_per_evf,
    replace = FALSE
)

gene_indices <- matrix(gene_indices, ncol = n_evfs_to_add)
# Add Components
# 	- Component 1:  Stratifies cells in half
# 	- Component 2:  Some effect of subset of memory cells
# 	- Component 3: Some effect of different subset of memory cells (maybe a small subset?)
#   - Component 4: Stratifies only the naïve cells in half
#   - Component 5: Stratifies all cells in half again

cell_indices <- list()
new_evfs <- list()

cell_indices[[1]] <- seq(ncells_total)
new_evfs[[1]] <- rnorm(ncells_total)*.5

mem_cells <- cell_indices[[1]][new_evfs[[1]] > 0]
naive_cells <- cell_indices[[1]][new_evfs[[1]] < 0]

cell_indices[[2]] <- sample(mem_cells, size=300)
new_evfs[[2]] <- rep(0, ncells_total)
new_evfs[[2]][cell_indices[[2]]] <- rnorm(
    length(cell_indices[[2]]))*.5 + 1

cell_indices[[3]] <- sample(mem_cells, size=30)
new_evfs[[3]] <- rep(0, ncells_total)
new_evfs[[3]][cell_indices[[3]]] <- rnorm(
    length(cell_indices[[3]]))*.5 + 1

cell_indices[[4]] <- naive_cells
new_evfs[[4]] <- rep(0, ncells_total)
new_evfs[[4]][cell_indices[[4]]] <- rnorm(
    length(cell_indices[[4]]))*.5 + 1

cell_indices[[5]] <- seq(ncells_total)
new_evfs[[5]] <- rnorm(ncells_total)*.2


param_names <- c("kon", "koff", "s")

for (i in seq(n_evfs_to_add)){
    evf_name <- paste0("new_evf", i)
    # Add the evf and initial gene_effect values
    cells_i <- cell_indices[[i]]
    for (p in seq(3)){
        evf_name_p <- paste0(param_names[p], "_", evf_name)
        new_vals <- new_evfs[[i]]
        evfs[[p]] <- cbind(evfs[[p]], new_vals)
        colnames(evfs[[p]])[ncol(evfs[[p]])] <- evf_name_p
        
        new_gvals <- rep(0, ngenes_total)
        gene_effects[[p]] <- cbind(gene_effects[[p]], new_gvals)
        colnames(gene_effects[[p]])[ncol(gene_effects[[p]])] <- evf_name_p
    }

    # Add the gene effects
    genes_i <- gene_indices[, i]

    # kon
    new_gvals <- rep(0, ngenes_total)
    new_gvals[genes_i] <- rnorm(length(genes_i))*new_effect_sd + new_effect_mean
    gene_effects[[1]][, ncol(gene_effects[[1]])] <- new_gvals

    # koff
    new_gvals <- rep(0, ngenes_total)
    new_gvals[genes_i] <- rnorm(length(genes_i))*new_effect_sd + new_effect_mean
    gene_effects[[2]][, ncol(gene_effects[[2]])] <- new_gvals

    # s
    new_gvals <- rep(0, ngenes_total)
    new_gvals[genes_i] <- rnorm(length(genes_i))*new_effect_sd + new_effect_mean
    gene_effects[[3]][, ncol(gene_effects[[3]])] <- new_gvals
}

# Simulate true counts
# true_counts_res has: counts, gene_effects, cell_meta, params
evf_res <- list(evfs = evfs, meta = meta)
true_counts_res <- trueCounts(evf_res, gene_effects)

true_counts <- true_counts_res[[1]]
cell_meta <- true_counts_res[[3]]
params <- true_counts_res[[4]]

# scaled_counts <- t(t(true_counts) / colSums(true_counts) * 10000)
# 
# tsne_res <- Rtsne(t(log1p(scaled_counts)), dims = 2, initial_dims = 30)

# ggplot() + aes(
#     x=tsne_res$Y[, 1],
#     y=tsne_res$Y[, 2],
#     #color=cell_meta$s_new_evf3
#     #color=evfs[[3]][, "s_new_evf1"]
#     color=seq(ncells_total) %in% cell_indices[, 1]
#     ) + geom_point()

# %% Why are some components just bad??
# Look at gene correlations

# g1 <- gene_indices[, 1]
# g11 <- g1[4]
# 
# gene_effects[[3]][g11, ]["s_new_evf1"]
# 
# plot(evfs[[3]][, "s_new_evf1"], scaled_counts[g11, ])
# 
# signal_genes <- log1p(scaled_counts)[as.numeric(gene_indices), ]
# gene_correlations <- cor(t(signal_genes))
# 
# X11.options(type="Xlib")
# 
# library(RColorBrewer)
# colors <- brewer.pal(11, name="RdBu")
# breaks <- seq(-.2, .2, length.out=12)
# heatmap(gene_correlations, Rowv=NA, Colv=NA, scale="none",
#     col=colors, breaks=breaks)


# What about the technical counts then

path <- system.file("data/gene_len_pool.RData", package = "SymSim")
load(path)
gene_len <- sample(gene_len_pool, ngenes_total, replace = FALSE)

options(mc.cores=20)
observed_counts_res <- True2ObservedCounts(
    true_counts = true_counts_res[[1]],
    meta_cell = true_counts_res[[3]],
    protocol = "UMI",
    alpha_mean = 0.05, alpha_sd = 0.002, lenslope = 0.01,
    nbins = 20, gene_len = gene_len,
    amp_bias_limit = c(-0.3, 0.3), rate_2PCR = 0.7,
    nPCR = 16, depth_mean = 5e4, depth_sd = 3e3
)

observed_counts <- observed_counts_res[[1]]
observed_meta_cell <- observed_counts_res[[2]]
nreads_perUMI <- observed_counts_res[[3]]
nUMI2seq <- observed_counts_res[[4]]

scaled_observed_counts <- t(t(observed_counts) / colSums(observed_counts) * 10000)

# Plot

# tsne_obs_res <- Rtsne(t(log1p(scaled_observed_counts)), dims = 2, initial_dims = 30)
# 
# ggplot() + aes(
#     x=tsne_obs_res$Y[, 1],
#     y=tsne_obs_res$Y[, 2],
#     color=observed_meta_cell$s_new_evf1
#     #color=colSums(observed_counts),
#     ) + geom_point()


# %% create output loom file

gene_ids <- paste0("Gene", seq(nrow(observed_counts)))
cell_ids <- paste0("Cell", seq(ncol(observed_counts)))

# Loom files broken in R, make it in python instead...
observed_counts_df <- data.frame(observed_counts)
rownames(observed_counts_df) <- gene_ids
colnames(observed_counts_df) <- cell_ids

write.table(
    observed_counts_df, gzfile(obs_counts_file),
    sep = "\t", col.names=NA)

# Other output files

# True counts

true_counts_df <- data.frame(true_counts)
rownames(true_counts_df) <- gene_ids
colnames(true_counts_df) <- cell_ids

write.table(true_counts_df, gzfile(true_counts_file), sep = "\t", col.names=NA)

# Component indices
gene_indices_df <- data.frame(gene_indices - 1) # -1 for 0-based
write.table(gene_indices_df, gene_indices_file, sep = "\t", col.names=NA)

# Component effects (for S)
gene_effects_df <- data.frame(gene_effects[[3]])
rownames(gene_effects_df) <- gene_ids
write.table(gene_effects_df, gene_effects_file, sep = "\t", col.names=NA)

cell_effects_df <- data.frame(evfs[[3]])
rownames(cell_effects_df) <- cell_ids
write.table(cell_effects_df, cell_effects_file, sep = "\t", col.names=NA)

