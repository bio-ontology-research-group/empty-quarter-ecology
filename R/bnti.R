#!/usr/bin/env Rscript
# β-NTI (beta-Nearest Taxon Index) and Raup-Crick following Stegen et al. 2012/2013.
#
# Arguments (positional):
#   1. FT_TSV   — ASV × sample feature table (TSV, first column is ASV id)
#   2. TREE_NWK — Newick tree covering the ASVs
#   3. META_TSV — sample metadata with columns 'sample' and 'compartment'
#   4. OUT_DIR  — writable output directory
#
# Environment:
#   N_NULL    — number of null iterations (default 999)
#   N_TOP     — top-prevalence ASVs per compartment (default 2000)
#   COMPARTMENTS — comma-separated list (default "surface,deep,rhizosphere")

suppressPackageStartupMessages({
  lib <- Sys.getenv("R_LIBS_USER", unset = NA)
  if (!is.na(lib) && nzchar(lib)) {
    dir.create(lib, recursive = TRUE, showWarnings = FALSE)
    .libPaths(c(lib, .libPaths()))
  }
  for (p in c("ape", "picante", "vegan", "data.table", "parallel")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org",
                       lib = .libPaths()[1], Ncpus = 4)
    }
  }
  library(ape); library(picante); library(vegan); library(data.table)
  library(parallel)
})

args     <- commandArgs(trailingOnly = TRUE)
ft_tsv   <- args[1]
tree_nwk <- args[2]
meta_tsv <- args[3]
out_dir  <- args[4]
n_null   <- as.integer(Sys.getenv("N_NULL", "999"))
n_top    <- as.integer(Sys.getenv("N_TOP",  "1000"))
n_cores  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
comps    <- strsplit(Sys.getenv("COMPARTMENTS",
                                "surface,deep,rhizosphere"), ",")[[1]]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()
log <- function(...) message(sprintf("[%s] ", Sys.time()), ...)

log("reading feature table ", ft_tsv)
ft <- fread(ft_tsv, sep = "\t", header = TRUE)
asv <- ft[[1]]; ft[[1]] <- NULL
ft <- as.matrix(ft); rownames(ft) <- asv
log(sprintf("  %d ASVs × %d samples", nrow(ft), ncol(ft)))

log("reading tree ", tree_nwk)
tree <- read.tree(tree_nwk)
log(sprintf("  %d tips", length(tree$tip.label)))

log("reading metadata ", meta_tsv)
meta <- fread(meta_tsv, sep = "\t")
stopifnot(all(c("sample", "compartment") %in% colnames(meta)))

for (comp in comps) {
  log("=== compartment: ", comp, " ===")
  samps <- meta[compartment == comp, sample]
  samps <- intersect(samps, colnames(ft))
  log(sprintf("  samples: %d", length(samps)))
  ft_c <- ft[, samps, drop = FALSE]
  prev <- rowSums(ft_c > 0)
  top_asvs <- names(sort(prev, decreasing = TRUE)[seq_len(min(n_top, length(prev)))])
  top_asvs <- intersect(top_asvs, tree$tip.label)
  ft_c <- ft_c[top_asvs, , drop = FALSE]
  log(sprintf("  ASVs kept: %d", nrow(ft_c)))

  tree_c <- keep.tip(tree, top_asvs)
  phy_d  <- cophenetic.phylo(tree_c)

  comm <- t(ft_c)
  comm <- comm / rowSums(comm)
  storage.mode(comm) <- "double"

  log("  observed β-MNTD (single-thread)")
  obs <- as.matrix(comdistnt(comm, phy_d, abundance.weighted = TRUE))
  log(sprintf("  observed done; %d null iterations on %d cores",
              n_null, n_cores))

  n_s <- nrow(comm)
  tip_labels_orig <- tree_c$tip.label
  sample_names <- rownames(comm)
  asv_names <- colnames(comm)

  # Pre-capture the objects each worker needs; mclapply inherits env on Linux
  run_one <- function(seed) {
    set.seed(seed)
    t2 <- tree_c
    t2$tip.label <- sample(tip_labels_orig)
    pd <- cophenetic.phylo(t2)[asv_names, asv_names]
    as.matrix(comdistnt(comm, pd, abundance.weighted = TRUE))
  }

  null_list <- mclapply(seq_len(n_null), run_one, mc.cores = n_cores,
                        mc.preschedule = FALSE)
  # mclapply returns error objects on worker failure — flag any
  ok <- vapply(null_list, is.matrix, logical(1))
  log(sprintf("  null iterations successful: %d / %d", sum(ok), n_null))
  stopifnot(all(ok))

  null_arr <- simplify2array(null_list)
  null_mean <- apply(null_arr, c(1, 2), mean)
  null_sd   <- apply(null_arr, c(1, 2), stats::sd)
  bnti      <- (obs - null_mean) / null_sd
  dimnames(bnti) <- dimnames(obs)
  dimnames(null_mean) <- dimnames(obs)
  dimnames(null_sd) <- dimnames(obs)

  fwrite(data.frame(obs),       file = file.path(out_dir, sprintf("obs_bMNTD_%s.tsv", comp)),
         sep = "\t", row.names = TRUE)
  fwrite(data.frame(null_mean), file = file.path(out_dir, sprintf("null_mean_%s.tsv", comp)),
         sep = "\t", row.names = TRUE)
  fwrite(data.frame(null_sd),   file = file.path(out_dir, sprintf("null_sd_%s.tsv",   comp)),
         sep = "\t", row.names = TRUE)
  fwrite(data.frame(bnti),      file = file.path(out_dir, sprintf("bNTI_%s.tsv",     comp)),
         sep = "\t", row.names = TRUE)

  v <- bnti[upper.tri(bnti)]
  log(sprintf("  β-NTI summary:  mean %.3f  sd %.3f  |β-NTI|>2 frac: %.3f",
              mean(v, na.rm = TRUE), sd(v, na.rm = TRUE),
              mean(abs(v) > 2, na.rm = TRUE)))
}

log("all compartments done")
