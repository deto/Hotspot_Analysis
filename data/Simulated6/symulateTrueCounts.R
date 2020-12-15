devtools::load_all("/data/yosef2/users/david.detomaso/repos/SymSim/")
library(Rtsne)
library(loomR)
library(ggplot2)
library(pbmcapply)

gene_effects_file <- snakemake@output[["gene_effects"]]
true_counts_file <- snakemake@output[["true_counts"]]
cell_meta_file <- snakemake@output[["cell_meta"]]

if ("CELL_FACTOR" %in% names(snakemake@params)) {
    CELL_FACTOR <- snakemake@params[["CELL_FACTOR"]]
} else {
    CELL_FACTOR <- 1
}

# Inside this file I've modified the SimulateTrueCounts
# Function to add an offset term to each gene

SimulateTrueCountsOffset <- function(ncells_total,min_popsize,i_minpop=1,ngenes, 
                               evf_center=1,evf_type="one.population",nevf=10,
                               n_de_evf=0,impulse=F,vary='all',
                               Sigma=0.5,phyla=NULL,geffect_mean=0,gene_effects_sd=1,gene_effect_prob=0.3,
                               bimod=0,param_realdata="zeisel.imputed",joint=F,randseed,SE=F){
  set.seed(randseed)
  n_nd_evf=nevf-n_de_evf
  seed <- sample(c(1:1e5),size=2)
  param_names <- c("kon", "koff", "s")
  if(evf_type=='one.population'){
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
  gene_effects <- GeneEffects(ngenes=ngenes,nevf=nevf,randseed=seed[2],prob=gene_effect_prob,
                              geffect_mean=geffect_mean,geffect_sd=gene_effects_sd)

  ############ Modification starts HERE
  evfs <- evf_res[[1]]
  meta <- evf_res[[2]]

  # Add an offset term
  evfs <- lapply(evfs, function(evf) {
      evf[, 1] <- 1.0
      return(evf)
  })

  # With random coupling to s, kon, and koff
  gene_effects <- lapply(gene_effects, function(ge) {
      ge[, 1] <- rnorm(nrow(ge), mean = geffect_mean, sd = gene_effects_sd)
      return(ge)
  })

  evf_res <- list(evfs, meta)
  ############## End Modification

  if(!is.null(param_realdata)){
    if(param_realdata=="zeisel.imputed"){
      path <- system.file('data/param_realdata.zeisel.imputed.robj', package = 'SymSim')
      load(path)
    }
    if(param_realdata=="zeisel.pop4"){
      path <- system.file('data/param_realdata.zeisel.pop4.robj', package = 'SymSim')
      load(path)
    }
    match_params[,1]=log(base=10,match_params[,1])
    match_params[,2]=log(base=10,match_params[,2])
    match_params_den <- lapply(c(1:3),function(i){
      density(match_params[,i],n=2000)
    })
    if(joint==F){
      params <- Get_params(gene_effects,evf_res[[1]],match_params_den,bimod)      
    }else{params <- Get_params3(gene_effects,evf_res[[1]],match_params,bimod)}
  }else{
    params <- Get_params2(gene_effects,evf_res[[1]],bimod,ranges)
  }
  
  counts <- lapply(c(1:ngenes),function(i){
    count <- sapply(c(1:ncells_total),function(j){
      y <- rbeta(1,params[[1]][i,j],params[[2]][i,j])
      x <- rpois(1,y*params[[3]][i,j])
      return(x)
    })
  })
  cell_meta <- cbind( cellid=paste('cell',seq(1,ncells_total),sep='_'),evf_res[[2]],evf_res[[1]])
  counts <- do.call(rbind,counts)
  return(list(counts,gene_effects,cell_meta,params))
}



# %% As in the vignette

phyla <- read.tree(
    text = "(((A:1,B:1):1,(C:0.5,D:0.5):1.5):1,E:3);"
)

true_counts_res <- SimulateTrueCountsOffset(
    ncells_total = 3000 * CELL_FACTOR,
    min_popsize = 200 * CELL_FACTOR,
    ngenes = 5000,
    evf_center = 1,
    evf_type = "discrete",
    nevf = 6,
    n_de_evf = 5,
    Sigma = 0.5,
    vary = "s",
    gene_effect_prob = 0.02,
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
