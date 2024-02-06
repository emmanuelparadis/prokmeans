## k_means.R (2024-02-06)

##   Probabilistic K-Means

## Copyright 2023-2024 Emmanuel Paradis

## This file is part of the R-package `prokmeans'.
## See the file ../COPYING for licensing issues.

pkm <- function(x, cls, K, threshold = 1e-5, iter.lim = 200, quiet = FALSE)
{
    if (is.vector(x)) x <- as.matrix(x)
    if (length(cls) != nrow(x))
        stop("number of values in 'cls' not equal to number of rows in 'x'")
    PARA <- as.integer(c(K, threshold * length(cls), iter.lim, quiet))
    .Call(C_kmeans_dnorm, x, cls, PARA)
}

skm <- function(x, cls, K, algorithm = "Hartigan-Wong", iter.lim = 10, quiet = FALSE)
{
    if (is.vector(x)) x <- as.matrix(x)
    if (missing(cls)) {
        if (missing(K)) stop("arguments 'cls' and 'K' cannot be both missing")
        if (length(K) > 1) {
            warning("argument 'K' is longer than one: taking its 1st value")
            K <- K[1]
        }
        cls <- K
    } else {
        if (length(cls) != nrow(x))
            stop("number of values in 'cls' not equal to number of rows in 'x'")
    }
    kmeans(x, cls, iter.lim, algorithm = algorithm, trace = quiet)$cluster
}

kmeansConstr <- function(x, cluster, K, iter.lim = 100,
                         restart.max = 10, diag4restart = 10,
                         quiet = FALSE)
{
    x <- as.matrix(x) # make sure x is a matrix
    if (nrow(x) != length(cluster))
        stop("nrow(x) and length(cluster) must be equal.")
    cls0 <- as.integer(cluster)
    cls0[is.na(cluster)] <- 0L
    cls0[cls0 < 0L] <- 0L
    mx <- max(cls0)
    if (mx > K)
        stop(sprintf("Maximum value in cluster (%d) greater than K (%d).\nMaybe you want K = %d.\n", mx, K, mx))
    PARS <- as.integer(c(K, iter.lim, restart.max, diag4restart, !quiet))
    .Call(kmeansConstr_Call, x, cls0, PARS)
}
