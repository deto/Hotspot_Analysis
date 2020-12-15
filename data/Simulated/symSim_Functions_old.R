devtools::load_all("/data/yosef2/users/david.detomaso/repos/SymSim/")

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

    counts <- lapply(c(1:ngenes), function(i){
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


generate_data <- function(ncells_total, ngenes_total,
                          n_genes_signal, signal_overlap_proportion,
                          nevf, min_popsize = 200, randseed = 1, i_minpop = 1,
                          evf_type = "one.population", n_de_evf = 0,
                          impulse = F, vary = "all", phyla = NULL,
                          Sigma_signal = 0.5
                          ){

    # Method
    n_genes_overlap <- floor(n_genes_signal * signal_overlap_proportion)

    geffect_mean <- 0
    geffect_sd <- 1
    evf_center <- 1
    Sigma <- 0.5

    rr <- compute_evf_params(
            ncells_total,
            ngenes = (ngenes_total - n_genes_signal + n_genes_overlap),
            geffect_mean = geffect_mean, gene_effects_sd = geffect_sd,
            evf_center = evf_center, Sigma = Sigma, nevf = nevf,
            min_popsize = min_popsize, randseed = randseed, i_minpop = i_minpop,
            evf_type = evf_type, n_de_evf = n_de_evf, impulse = impulse,
            vary = vary, phyla = phyla
    )

    evf_res <- rr[[1]]
    evfs <- evf_res[[1]]
    meta <- evf_res[[2]]

    gene_effects <- rr[[2]]

    s <- matrix(rnorm(ncells_total, evf_center, Sigma_signal), ncol = 1)
    zeros <- s * 0
    evf_signal <- list(zeros, zeros, s)


    gene_effects_signal <- GeneEffects(
        ngenes = n_genes_signal, nevf = 1, randseed = 1000, prob = 1,
        geffect_mean = geffect_mean, geffect_sd = geffect_sd
    )

    prepend_N_genes <- function(ge, N) {
        result <- lapply(ge, function(x){
            y <- matrix(0, nrow = N, ncol = ncol(x))
            z <- rbind(y, x)
            return(z)
        })

        return(result)
    }

    append_N_genes <- function(ge, N) {
        result <- lapply(ge, function(x){
            y <- matrix(0, nrow = N, ncol = ncol(x))
            z <- rbind(x, y)
            return(z)
        })

        return(result)
    }

    gene_effects <- prepend_N_genes(gene_effects, n_genes_signal - n_genes_overlap)
    gene_effects_signal <- append_N_genes(gene_effects_signal, ngenes_total - n_genes_signal)

    # combine evfs
    combine_evfs <- function(evfA, evfB){
        evf_all <- lapply(1:3, function(i){
            evfAi <- evfA[[i]]
            evfBi <- evfB[[i]]
            evfi <- cbind(evfAi, evfBi)
            return(evfi)
        })
    }

    combine_gene_effects <- function(geA, geB){
        ge_all <- lapply(1:3, function(i){
            geAi <- geA[[i]]
            geBi <- geB[[i]]
            gei <- cbind(geAi, geBi)
            return(gei)
        })
    }

    evf_all <- combine_evfs(evfs, evf_signal)
    gene_effects_all <- combine_gene_effects(gene_effects, gene_effects_signal)


    evf_res <- list(evfs = evf_all, meta = meta)
    gene_effects <- gene_effects_all

    # out is a list(counts,gene_effects,cell_meta,params)
    true_counts_res <- trueCounts(evf_res, gene_effects)

    path <- system.file("data/gene_len_pool.RData", package = "SymSim")
    load(path)
    gene_len <- sample(gene_len_pool, ngenes_total, replace = FALSE)


    observed_counts <- True2ObservedCounts(
        true_counts = true_counts_res[[1]],
        meta_cell = true_counts_res[[3]],
        protocol = "UMI",
        alpha_mean = 0.05, alpha_sd = 0.02, lenslope = 0.01,
        nbins = 20, gene_len = gene_len,
        amp_bias_limit = c(-0.3, 0.3), rate_2PCR = 0.7,
        nPCR = 16, depth_mean = 5e4, depth_sd = 3e3)


    meta_cell <- observed_counts[[2]]
    counts <- observed_counts[[1]]

    return(list(counts, evf_res[[1]], gene_effects, meta_cell))
}

# Fixed Parameters
ncells_total <- 4000
ngenes_total <- 500

# Variable Parameters

params_file <- snakemake@input[["params"]]
job_number <- snakemake@params[["job_number"]]

params <- read.table(params_file, sep = "\t", header = TRUE)
params <- params[params$job == job_number, ]

n_genes_signal <- params[, "n_genes_signal"] # Number of genes in the signal (default 100)
signal_overlap_proportion <- params[, "signal_overlap"] # proportion of overlap between DE genes and signal genes (default 0.0)
signal_strength <- params[, "signal_strength"]  # std-deviation of evf for signal, 0.5 is default

message("Generating counts with parameters...")
message("n_genes_signal: ", n_genes_signal)
message("signal_overlap: ", signal_overlap_proportion)
message("signal_strength: ", signal_strength)

# Outputs
out_counts <- snakemake@output[["counts"]]
out_evf <- snakemake@output[["signal_evf"]]
out_meta <- snakemake@output[["meta"]]

out_dir <- dirname(out_counts)
dir.create(out_dir, showWarnings = FALSE)

# ok, so given evf_res and gene_effects we need to augment these with our signal
nevf <- 9

min_popsize <- 2000
evf_type <- "discrete"
n_de_evf <- 8

tree <- "(A:5,B:5);"
phyla <- read.tree(text = tree)


rr <- generate_data(ncells_total, ngenes_total,
                    n_genes_signal,
                    signal_overlap_proportion,
                    nevf, min_popsize = min_popsize,
                    evf_type = evf_type, n_de_evf = n_de_evf,
                    phyla = phyla, vary = "s", Sigma_signal = signal_strength)

counts <- rr[[1]]
evfs <- rr[[2]]
gene_effects <- rr[[3]]
meta_cell <- rr[[4]]

# Save results to file
library(feather)
counts <- as.data.frame(counts)
write_feather(counts, out_counts)

write.table(evfs[[3]][, nevf + 1, drop = F], out_evf, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table(meta_cell[, c("cellid", "pop", "depth")], out_meta, sep = "\t",
            row.names = FALSE, col.names = TRUE)
