#!/usr/bin/env Rscript
# Generate fresh 100-Mb chromosomes using scrm's sequential coalescent approximation.
# One chromosome at a time limits peak memory. Requires CRAN scrm 1.7.5.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L, requireNamespace("scrm", quietly=TRUE))
out <- args[1]
dir.create(out,showWarnings=FALSE,recursive=TRUE)
scenarios <- data.frame(scenario=1:4, Na=c(18000,24000,15000,28000),
                        fraction=c(.68,.46,.60,.80), T=c(1200,3200,6500,2800))
write.table(scenarios,file.path(out,"scenarios.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
L <- 100000000
n <- 114L
reps <- 10L
# The 10-kb pruning window controls the sequential approximation, not the
# physical chromosome length. No independent segment resets are introduced.
simulate_one <- function(s, rep) {
  p <- scenarios[s,]
  current <- p$Na*p$fraction
  set.seed(465201L+100*s+rep)
  opts <- sprintf("%d 1 -r %.12g %d -t %.12g -l 10000 -SC abs -eN %.12g %.12g",
                  n,4*current*1.87e-8*(L-1),L,4*current*3.42e-8*L,
                  p$T/(4*current),1/p$fraction)
  mat <- scrm::scrm(opts)$seg_sites[[1]]
  pos <- as.numeric(colnames(mat))
  ac <- as.integer(colSums(mat))
  stopifnot(nrow(mat)==n,all(is.finite(pos)),all(pos>=0 & pos<L),
            all(ac>=1L & ac<n))
  block <- floor(pos/1000000)
  x <- matrix(tabulate(block*(n-1L)+ac, nbins=100L*(n-1L)),nrow=100L,byrow=TRUE)
  stopifnot(sum(x)==length(pos))
  target <- file.path(out,sprintf("sim_%d_%02d.tsv",s,rep))
  stage <- paste0(target,".partial")
  write.table(x,stage,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  stopifnot(file.rename(stage,target))
}
for (r in seq_len(reps)) for (s in scenarios$scenario) {
  simulate_one(s,r)
  gc(verbose=FALSE)
}
