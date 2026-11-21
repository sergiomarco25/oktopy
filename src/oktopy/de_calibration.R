#!/usr/bin/env Rscript
# Shuffled-label calibration of the DE cache (oktopy engine, phase 7a): for a random sample of the dataset's
# region x class x contrast fits, shuffle the case / reference labels among the fit's units (covariates stay with their
# units) and refit with the same model as de_cache.R. Under shuffled labels no gene is truly different, so a calibrated
# test calls about 5% of genes at p < 0.05 and rarely any at FDR < 0.05.
# Usage: Rscript de_calibration.R <dataset_dir> <out_dir> <perms> <max_fits> <seed>
#   dataset_dir = results/de_cache/<dataset> (input/manifest.tsv, units.tsv, counts); writes <out_dir>/de_null.tsv and
#   <out_dir>/null/<region>__<class>__<contrast>__p<k>.tsv.gz (gene, logFC, p) for the gene-set calibration.
# fit_one is a verbatim copy of de_cache.R's (tests/test_robustness.py fails if they differ).
suppressPackageStartupMessages({ library(data.table); library(edgeR) })
args <- commandArgs(trailingOnly = TRUE)
dsdir <- args[1]; outdir <- args[2]
perms <- as.integer(args[3]); max_fits <- as.integer(args[4]); seed <- as.integer(args[5])
indir <- file.path(dsdir, "input"); dir.create(file.path(outdir, "null"), recursive = TRUE, showWarnings = FALSE)
units <- fread(file.path(indir, "units.tsv"), colClasses = "character")
man <- fread(file.path(indir, "manifest.tsv"))
scl <- function(v) as.numeric(scale(v))

fit_one <- function(counts, u, term) {
  covs <- c()
  if (all(!is.na(suppressWarnings(as.numeric(u$age)))) && sd(as.numeric(u$age)) > 0) { u$age_s <- scl(as.numeric(u$age)); covs <- c(covs, "age_s") }
  if (all(u$sex %in% c("M", "F")) && length(unique(u$sex)) > 1) { u$sex <- factor(u$sex); covs <- c(covs, "sex") }
  if (all(!is.na(suppressWarnings(as.numeric(u$pmi)))) && sd(as.numeric(u$pmi)) > 0) { u$pmi_s <- scl(as.numeric(u$pmi)); covs <- c(covs, "pmi_s") }
  # keep covariates only while the design stays estimable with >= 2 residual df
  repeat {
    design <- model.matrix(as.formula(paste("~", paste(c(term, covs), collapse = " + "))), data = u)
    if (qr(design)$rank == ncol(design) && nrow(design) - ncol(design) >= 2) break
    if (!length(covs)) return(NULL)
    covs <- head(covs, -1)
  }
  y <- DGEList(counts); keep <- filterByExpr(y, design)
  if (sum(keep) < 50) return(NULL)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design); fit <- glmQLFit(y, design, robust = TRUE)
  coef <- grep(paste0("^", term), colnames(design), value = TRUE)[1]
  tt <- as.data.table(topTags(glmQLFTest(fit, coef = coef), n = Inf)$table, keep.rownames = "gene")
  setnames(tt, c("PValue", "FDR"), c("p", "padj"))
  tt[, covariates := paste(covs, collapse = "+")]
  tt
}

set.seed(seed)
todo <- unique(man[run == TRUE, .(region, cell_class, file, contrast)])
if (nrow(todo) > max_fits) todo <- todo[sort(sample.int(nrow(todo), max_fits))]
rows <- list()
for (fl in unique(todo$file)) {
  counts <- fread(cmd = paste("gzip -dc", shQuote(file.path(indir, fl))))
  gn <- counts[[1]]; counts <- as.matrix(counts[, -1, with = FALSE]); rownames(counts) <- gn
  u <- units[match(colnames(counts), units$unit)]
  stem <- sub("\\.csv\\.gz$", "", fl)
  for (j in which(todo$file == fl)) {
    con <- todo$contrast[j]
    g <- u[[con]]; sel <- g %in% c("case", "ref")
    for (k in seq_len(perms)) {
      uu <- u[sel]; uu$group <- factor(sample(g[sel]), levels = c("ref", "case"))
      tt <- tryCatch(fit_one(counts[, sel, drop = FALSE], uu, "group"), error = function(e) NULL)
      if (is.null(tt)) next
      fwrite(tt[, .(gene, logFC, p)], file.path(outdir, "null", paste0(stem, "__", con, "__p", k, ".tsv.gz")), sep = "\t")
      chi <- qchisq(pmax(tt$p, 1e-300), 1, lower.tail = FALSE)
      rows[[length(rows) + 1]] <- data.table(file = fl, region = todo$region[j], cell_class = todo$cell_class[j],
        contrast = con, perm = k, n_case = sum(uu$group == "case"), n_ref = sum(uu$group == "ref"),
        n_genes = nrow(tt), frac_p05 = mean(tt$p < 0.05), n_fdr05 = sum(tt$padj < 0.05),
        lambda = median(chi) / qchisq(0.5, 1), covariates = tt$covariates[1])
    }
  }
}
res <- rbindlist(rows)
fwrite(res, file.path(outdir, "de_null.tsv"), sep = "\t")
cat(nrow(res), "null fits;", "mean share of genes at p < 0.05:", round(mean(res$frac_p05), 4), "\n")
