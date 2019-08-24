library(idr)

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

    IDR.level <- 1e-2
    x.selected <- select.IDR(x, res$IDR, IDR.level)

    x.selected$n

    res <- c(fileA, fileB, x.selected$n)
    return(res)
}

files <- c(
    "Puck_180819_9/hotspot/hotspot.txt",
    "Puck_180819_10/hotspot/hotspot.txt",
    "Puck_180819_11/hotspot/hotspot.txt",
    "Puck_180819_12/hotspot/hotspot.txt"
)

stat_col <- "Z"

combos <- make_pairs(files)

results <- lapply(combos, pair_fun, stat_col)
results <- as.data.frame(do.call(rbind, results))
colnames(results) <- c("fileA", "fileB", "n")
results_hs <- results

files <- c(
    "Puck_180819_9/spatialDE/spatialDE.txt",
    "Puck_180819_10/spatialDE/spatialDE.txt",
    "Puck_180819_11/spatialDE/spatialDE.txt",
    "Puck_180819_12/spatialDE/spatialDE.txt"
)

stat_col <- "LLR"

combos <- make_pairs(files)

results <- lapply(combos, pair_fun, stat_col)
results <- as.data.frame(do.call(rbind, results))
colnames(results) <- c("fileA", "fileB", "n")
results_sde <- results


combos <- list(
    c("Puck_180819_9/hotspot/hotspot.txt",  "Puck_180819_9/spatialDE/spatialDE.txt"),
    c("Puck_180819_10/hotspot/hotspot.txt", "Puck_180819_10/spatialDE/spatialDE.txt"),
    c("Puck_180819_11/hotspot/hotspot.txt", "Puck_180819_11/spatialDE/spatialDE.txt"),
    c("Puck_180819_12/hotspot/hotspot.txt", "Puck_180819_12/spatialDE/spatialDE.txt")
)

stat_colA <- "Z"
stat_colB <- "LLR"

results <- lapply(combos, pair_fun, stat_colA, stat_colB)
results <- as.data.frame(do.call(rbind, results))
colnames(results) <- c("fileA", "fileB", "n")
results_comp <- results

# Some plots

fileA <- "Puck_180819_11/hotspot/hotspot.txt"
fileB <- "Puck_180819_12/hotspot/hotspot.txt"
stat_col <- "Z"

dfA <- read.table(fileA, sep = "\t", header = TRUE, row.names = 1)
dfB <- read.table(fileB, sep = "\t", header = TRUE, row.names = 1)

common <- intersect(row.names(dfA), row.names(dfB))

valsA <- dfA[common, stat_col]
valsB <- dfB[common, stat_col]

x <- cbind(valsA, valsB)

mu <- 1
sigma <- 0.5
rho <- 0.5
p <- 0.1
res <- est.IDR(x, mu, sigma, rho, p, eps = 0.001, max.ite = 500)

IDR.level <- 1e-2
x.selected <- select.IDR(x, res$IDR, IDR.level)

plot(x[, 1], x[, 2],
    xlim=c(-5, 100), ylim=c(-5, 100), cex=.2)
points(x.selected$x[, 1], x.selected$x[, 2], col = "red", cex=.2)


fileA <- "Puck_180819_9/spatialDE/spatialDE.txt"
fileB <- "Puck_180819_12/spatialDE/spatialDE.txt"
stat_col <- "LLR"

dfA <- read.table(fileA, sep = "\t", header = TRUE, row.names = 1)
dfB <- read.table(fileB, sep = "\t", header = TRUE, row.names = 1)

common <- intersect(row.names(dfA), row.names(dfB))

valsA <- dfA[common, stat_col]
valsB <- dfB[common, stat_col]

x <- cbind(valsA, valsB)

mu <- 1
sigma <- 0.5
rho <- 0.5
p <- 0.1
res <- est.IDR(x, mu, sigma, rho, p, eps = 0.001, max.ite = 500)

IDR.level <- 1e-2
x.selected <- select.IDR(x, res$IDR, IDR.level)

plot(x[, 1], x[, 2],
    xlim=c(-5, 100), ylim=c(-5, 100), cex=.2)
points(x.selected$x[, 1], x.selected$x[, 2], col = "red", cex=.2)

# What if we IDR all samples together?

files <- c(
    "Puck_180819_9/hotspot/hotspot.txt",
    "Puck_180819_10/hotspot/hotspot.txt",
    "Puck_180819_11/hotspot/hotspot.txt",
    "Puck_180819_12/hotspot/hotspot.txt"
)

stat_col <- "Z"
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

IDR.level <- 1e-2
x.selected <- select.IDR(x, res$IDR, IDR.level)

x.selected$n

files <- c(
    "Puck_180819_9/spatialDE/spatialDE.txt",
    "Puck_180819_10/spatialDE/spatialDE.txt",
    "Puck_180819_11/spatialDE/spatialDE.txt",
    "Puck_180819_12/spatialDE/spatialDE.txt"
)

stat_col <- "LLR"
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

IDR.level <- 1e-2
x.selected <- select.IDR(x, res$IDR, IDR.level)

x.selected$n
