#!/usr/bin/env Rscript
# Usage: Rscript --vanilla scripts/simulate.R DIRECTORY
out <- commandArgs(trailingOnly = TRUE)
stopifnot(length(out) == 1L, requireNamespace("phyclust", quietly = TRUE))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
scenarios <- data.frame(scenario=1:4, Na=c(20000, 20000, 25000, 15000),
                        fraction=c(0.75, 0.50, 0.30, 0.70), T=c(1600, 5000, 2500, 8000))
n <- 114L
mu <- 3.42e-8
recomb <- 1.87e-8
locus_bp <- 25000L
n_loci <- 1000L
n_blocks <- 100L
replicates <- 20L
write.table(scenarios, file.path(out, "scenarios.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

one_genome <- function(id) {
  scenario <- 1L + (id-1L) %/% replicates
  sample_id <- 1L + (id-1L) %% replicates
  model <- scenarios[scenario,]
  Nc <- model$Na * model$fraction
  opts <- sprintf("-t %.12g -r %.12g %d -eN %.12g %.12g",
                  4*Nc*mu*locus_bp, 4*Nc*recomb*(locus_bp-1), locus_bp,
                  model$T/(4*Nc), 1/model$fraction)
  # phyclust uses R's RNG. Each whole 25 Mb replicate receives a distinct seed.
  set.seed(17471L + id)
  counts <- matrix(0L, n_blocks, n-1L)
  batch <- 50L
  for (first in seq.int(1L, n_loci, batch)) {
    x <- phyclust::ms(n, batch, opts)
    starts <- which(startsWith(x, "segsites:"))
    stopifnot(length(starts) == batch)
    for (k in seq_len(batch)) {
      h <- starts[k]
      sites <- as.integer(sub("segsites: *", "", x[h]))
      if (!sites) next
      hap <- x[h + 1L + seq_len(n)]
      bits <- as.integer(charToRaw(paste0(hap, collapse=""))) - 48L
      stopifnot(length(bits) == n * sites)
      ac <- rowSums(matrix(bits, nrow=sites))
      block <- 1L + (first+k-2L) %/% (n_loci %/% n_blocks)
      counts[block,] <- counts[block,] + tabulate(ac, nbins=n-1L)
    }
  }
  path <- file.path(out, sprintf("sim_%d_%02d.tsv", scenario, sample_id))
  write.table(counts, path, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  TRUE
}
jobs <- parallel::mclapply(seq_len(nrow(scenarios)*replicates), one_genome,
                           mc.cores=3L, mc.set.seed=FALSE)
stopifnot(all(vapply(jobs, identical, logical(1), TRUE)))
