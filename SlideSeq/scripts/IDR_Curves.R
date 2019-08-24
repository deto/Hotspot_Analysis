library(idr)
library(ggplot2)
library(jsonlite)

make_pairs <- function(items){
    combos <- list()
    for (i in seq_along(items)){
        for (j in seq_along(items)){
            if (i >= j){
                next
            }
            combos[[length(combos) + 1]] <- c(items[i], items[j])
        }
    }
    return(combos)
}

pair_fun <- function(file_pair, stat_col, stat_colb=stat_col){
    fileA <- file_pair[1]
    fileB <- file_pair[2]

    dfA <- read.table(fileA, sep = "\t", header = TRUE, row.names = 1)
    dfB <- read.table(fileB, sep = "\t", header = TRUE, row.names = 1)

    common <- intersect(row.names(dfA), row.names(dfB))

    valsA <- dfA[common, stat_col]
    valsB <- dfB[common, stat_colb]

    x <- cbind(valsA, valsB)

    mu <- 1
    sigma <- 0.5
    rho <- 0.5
    p <- 0.1
    res <- est.IDR(x, mu, sigma, rho, p, eps = 0.001, max.ite = 500)
    res$fileA <- fileA
    res$fileB <- fileB

    return(res)
}

# %% Run for Hotspot pairs

files <- c(
    "../Puck_180819_9/hotspot/hotspot.txt",
    "../Puck_180819_10/hotspot/hotspot.txt",
    "../Puck_180819_11/hotspot/hotspot.txt",
    "../Puck_180819_12/hotspot/hotspot.txt"
)

stat_col <- "Z"
combos <- make_pairs(files)
results <- lapply(combos, pair_fun, stat_col)
jsonlite::write_json(
    results, "IDR_hotspot.json", auto_unbox = TRUE, pretty = TRUE, digits = NA)


# %% Run for SDE Pairs

files <- c(
    "../Puck_180819_9/spatialDE/spatialDE.txt",
    "../Puck_180819_10/spatialDE/spatialDE.txt",
    "../Puck_180819_11/spatialDE/spatialDE.txt",
    "../Puck_180819_12/spatialDE/spatialDE.txt"
)

stat_col <- "LLR"

combos <- make_pairs(files)
results <- lapply(combos, pair_fun, stat_col)
jsonlite::write_json(
    results, "IDR_SDE.json", auto_unbox = TRUE, pretty = TRUE, digits = NA)

# %% Pairs


combos <- list(
    c("../Puck_180819_9/hotspot/hotspot.txt",  "../Puck_180819_9/spatialDE/spatialDE.txt"),
    c("../Puck_180819_10/hotspot/hotspot.txt", "../Puck_180819_10/spatialDE/spatialDE.txt"),
    c("../Puck_180819_11/hotspot/hotspot.txt", "../Puck_180819_11/spatialDE/spatialDE.txt"),
    c("../Puck_180819_12/hotspot/hotspot.txt", "../Puck_180819_12/spatialDE/spatialDE.txt")
)

stat_colA <- "Z"
stat_colB <- "LLR"

results <- lapply(combos, pair_fun, stat_colA, stat_colB)
jsonlite::write_json(
    results, "IDR_Pairs.json", auto_unbox = TRUE, pretty = TRUE, digits = NA)

# %% Some plots

# What if we IDR all samples together?

idr_group <- function(files, stat_col) {

    dfs <- lapply(files, function(fileA) {
        read.table(fileA, sep = "\t", header = TRUE, row.names = 1)
    })

    common <- Reduce(
        intersect, lapply(dfs, function(x) row.names(x)),
        init = row.names(dfs[[1]])
    )

    vals <- lapply(dfs, function(df){
        df[common, stat_col]
    })

    x <- do.call(cbind, vals)

    mu <- 1
    sigma <- 0.5
    rho <- 0.5
    p <- 0.1
    res <- est.IDR(x, mu, sigma, rho, p, eps = 0.001, max.ite = 500)

    return(res)

}

files <- c(
    "../Puck_180819_9/hotspot/hotspot.txt",
    "../Puck_180819_10/hotspot/hotspot.txt",
    "../Puck_180819_11/hotspot/hotspot.txt",
    "../Puck_180819_12/hotspot/hotspot.txt"
)

stat_col <- "Z"

res_hs4 <- idr_group(files, stat_col)
res_hs3 <- idr_group(files[2:4], stat_col)

files <- c(
    "../Puck_180819_9/spatialDE/spatialDE.txt",
    "../Puck_180819_10/spatialDE/spatialDE.txt",
    "../Puck_180819_11/spatialDE/spatialDE.txt",
    "../Puck_180819_12/spatialDE/spatialDE.txt"
)

stat_col <- "LLR"
res_sd4 <- idr_group(files, stat_col)
res_sd3 <- idr_group(files[2:4], stat_col)


results <- list(
    "hs4" =  res_hs4,
    "hs3" =  res_hs3,
    "sd4" =  res_sd4,
    "sd3" =  res_sd3
)

jsonlite::write_json(
    results, "IDR_multi.json", auto_unbox = TRUE, pretty = TRUE, digits = NA)
