#!/usr/bin/env Rscript
# Fit a diploid, one-population instantaneous decline to a one-dimensional dadi SFS.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 4L) || !args[3] %in% c("folded", "polarized")) {
  stop("Usage: Rscript fit.R INPUT.fs OUTPUT_DIRECTORY folded|polarized [callable_length_bp]")
}
input <- args[1]
out <- args[2]
mode <- args[3]
mu <- 3.42e-8
region_bp <- if (length(args) == 4L) as.numeric(args[4]) else 100000000
if (!is.finite(region_bp) || region_bp <= 0) stop("Invalid callable length.")
lines <- readLines(input, warn = FALSE)
lines <- lines[nzchar(trimws(lines)) & !startsWith(trimws(lines), "#")]
header <- strsplit(trimws(lines[1]), "[[:space:]]+")[[1]]
n <- as.integer(header[1]) - 1L
if (length(lines) != 3L || n < 5L || !identical(header[2], "unfolded")) {
  stop("Expected a single-population unfolded dadi file, including its mask.")
}
obs <- scan(text = lines[2], quiet = TRUE)
mask <- scan(text = lines[3], quiet = TRUE)
i <- seq_len(n - 1L)
if (length(obs) != n + 1L || length(mask) != length(obs) ||
    any(!is.finite(obs)) || any(obs < 0) || any(mask[i + 1L] != 0)) {
  stop("Incorrect dimensions, invalid counts or unsupported masked interior bins.")
}
collapse_bins <- function(x) {
  if (mode == "polarized") return(x)
  # Each minor allele count has two possible derived counts, except n/2.
  as.numeric(rowsum(x, pmin(i, n - i), reorder = FALSE))
}
counts <- collapse_bins(obs[i + 1L])
total <- sum(counts)
if (total <= 0) stop("The spectrum contains no segregating sites.")

# In units of 2*N_current generations the neutral SFS generator is L %*% D.
# D_ii = i*(n-i)/2 and L has diagonal -2 and neighboring entries 1.
# Similarity transformation to sqrt(D) %*% L %*% sqrt(D) makes it symmetric.
drift <- i * (n - i) / 2
root <- sqrt(drift)
G <- diag(-2 * drift)
adjacent <- cbind(seq_len(n - 2L), 2:(n - 1L))
G[adjacent] <- G[adjacent[, 2:1, drop = FALSE]] <- root[-1L] * root[-length(root)]
eig <- eigen(G, symmetric = TRUE)
eq <- 1 / i
coeff <- drop(crossprod(eig$vectors, root * eq))
equilibrium <- collapse_bins(eq)
old_spectrum <- function(log_tau) {
  a <- exp(eig$values * exp(log_tau)) * coeff
  collapse_bins(pmax(0, drop(eig$vectors %*% a) / root))
}
# Profile the ancestral theta from the Poisson composite likelihood. Conditional
# on the total S it depends only on the shape, which simplifies the search.
positive <- counts > 0
deviance <- function(expected) {
  2 * sum(counts[positive] * log(counts[positive] / expected[positive]))
}
log_fraction_range <- log(c(0.0001, 1))
log_tau_range <- log(c(1e-5, 20))
profile <- function(log_tau, details = FALSE) {
  old <- old_spectrum(log_tau)
  evaluate <- function(log_fraction, details = FALSE) {
    fraction <- exp(log_fraction)
    shape <- fraction * equilibrium + (1 - fraction) * old
    theta <- total / sum(shape)
    expected <- theta * shape
    if (!details) return(deviance(expected))
    ancestral <- theta / (4 * mu * region_bp)
    current <- ancestral * fraction
    list(deviance = deviance(expected), ancestral_Ne = ancestral, current_Ne = current,
         generations = 2 * current * exp(log_tau), post_decline_fraction = fraction,
         log_tau = log_tau, expected = expected)
  }
  best <- optimize(evaluate, log_fraction_range, tol = 1e-9)
  candidates <- c(best$minimum, log_fraction_range)
  evaluate(candidates[which.min(vapply(candidates, evaluate, numeric(1)))], details)
}
grid <- seq(log_tau_range[1], log_tau_range[2], length.out = 65)
grid_values <- vapply(grid, profile, numeric(1))
idx <- which.min(grid_values)
interval <- grid[c(max(idx - 1L, 1L), min(idx + 1L, length(grid)))]
opt <- optimize(profile, interval, tol = 1e-9)
points <- c(grid[idx], opt$minimum)
fit <- profile(points[which.min(vapply(points, profile, numeric(1)))], TRUE)
if (idx %in% c(1, length(grid)) || fit$post_decline_fraction < 0.00011 ||
    fit$post_decline_fraction > 0.9999) {
  stop("Fit reached a search boundary or a constant-size solution; change time is unresolved.")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(count = seq_along(counts), observed = counts, expected = fit$expected),
            file.path(out, "spectrum.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
metrics <- data.frame(mode = mode, haploid_n = n, segregating_sites = total,
                      ancestral_Ne = fit$ancestral_Ne, current_Ne = fit$current_Ne,
                      generations = fit$generations, post_decline_fraction = fit$post_decline_fraction,
                      model_deviance = fit$deviance,
                      equilibrium_deviance = deviance(total * equilibrium / sum(equilibrium)))
write.table(metrics, file.path(out, "fit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(metrics, digits = 7, row.names = FALSE)
