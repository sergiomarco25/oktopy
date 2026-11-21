#!/usr/bin/env Rscript
# Sensitivity DE methods of the oktopy engine (phase 7b, method multiverse). Each runs on the DE cache's inputs
# (results/de_cache/<dataset>/input), with the DE cache's design and gene filter, so only the test changes:
#   edger_ql    the primary's model again (edgeR glmQLFit, robust); fit reproduces the DE cache, null adds the stage
#   limma_voom  voom precision weights, lmFit, eBayes(robust = TRUE)          ~ term + age + sex + pmi
#   deseq2      DESeq2 Wald test on the same model matrix, no independent filtering (filterByExpr already filtered)
#   wilcoxon    Wilcoxon rank-sum per gene on TMM-normalised CPM (Li Y 2022); a continuous stage: Spearman; no
#               covariates; logFC = difference of mean log-CPM (Spearman: rho)
# limma_voom with a sets file (oktopy.multiverse sets) also runs camera on the claims' gene sets of each fit: genes
# the claim expects down are sign-flipped first, so "Up" = as claimed; inter-gene correlation estimated per set (p, as
# Wu & Smyth 2012) and preset to 0.01 (p_preset, limma's default since 3.29.6).
# Usage:
#   Rscript de_sensitivity.R METHOD DATASET_DIR OUT_DIR SETS|- fit 0 0 0 [STAGE_CONTRAST STAGE_COLUMNS GROUP_COLUMN ZERO_GROUP]
#   Rscript de_sensitivity.R METHOD DATASET_DIR OUT_DIR SETS|- null PERMS MAX_FITS SEED [the same stage settings]
# fit   OUT_DIR/de_<method>/<region>__<class>__<contrast>.tsv.gz (gene, logFC, p, padj, stat, covariates),
#       OUT_DIR/de_log_<method>.tsv, and for limma_voom OUT_DIR/camera.tsv
# null  case / reference labels (or stage values) shuffled among a fit's units, covariates staying with their units:
#       OUT_DIR/null/de_null_<method>.tsv (one row per null fit), OUT_DIR/null/fdr_hits_<method>.tsv (genes at
#       FDR < 0.05 under the null), and for limma_voom OUT_DIR/null/camera_null.tsv
# design_for keeps the covariate rule of de_cache.R's fit_one verbatim (tests/test_robustness.py fails if they differ).
suppressPackageStartupMessages({ library(data.table); library(edgeR) })
args <- commandArgs(trailingOnly = TRUE)
method <- args[1]; dsdir <- args[2]; outdir <- args[3]; sets_file <- args[4]; mode <- args[5]
perms <- as.integer(args[6]); max_fits <- as.integer(args[7]); seed <- as.integer(args[8])
stage_contrast <- if (length(args) >= 10) args[9] else ""
stage_cols <- if (length(args) >= 10) strsplit(args[10], ",", fixed = TRUE)[[1]] else character()
group_col <- if (length(args) >= 11) args[11] else ""
zero_group <- if (length(args) >= 12) args[12] else ""
stopifnot(method %in% c("edger_ql", "limma_voom", "deseq2", "wilcoxon"), mode %in% c("fit", "null"))
ncpu <- max(1L, suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))), na.rm = TRUE)
if (method == "deseq2") {  # DESeq2 gives the same results in parallel
  suppressPackageStartupMessages(library(DESeq2))
  if (ncpu > 1) BiocParallel::register(BiocParallel::MulticoreParam(ncpu))
}
indir <- file.path(dsdir, "input")
units <- fread(file.path(indir, "units.tsv"), colClasses = "character")
man <- fread(file.path(indir, "manifest.tsv"))
sets <- if (sets_file != "-" && file.exists(sets_file) && file.size(sets_file) > 0)
  fread(sets_file, colClasses = list(character = c("stem", "contrast", "set_key", "gene"), numeric = "sign")) else NULL
scl <- function(v) as.numeric(scale(v))
stage_of <- function(u) {  # as oktopy.contrasts.stage_value
  b <- rep(NA_real_, nrow(u))
  for (col in stage_cols) if (col %in% names(u)) {
    v <- suppressWarnings(as.numeric(u[[col]]))
    if (any(!is.na(v))) { b <- v; break }
  }
  if (any(!is.na(b)) && nzchar(zero_group) && group_col %in% names(u)) b[is.na(b) & u[[group_col]] == zero_group] <- 0
  b
}

design_for <- function(u, term) {
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
  list(design = design, covs = covs, u = u)
}

camera_sets <- function(v, design, coef, s) {
  ci <- match(coef, colnames(design)); out <- list()
  for (k in unique(s$set_key)) {
    ss <- s[set_key == k & gene %in% rownames(v$E)]
    if (nrow(ss) < 2) next
    idx <- match(ss$gene, rownames(v$E)); flip <- idx[ss$sign < 0]
    if (length(flip)) v$E[flip, ] <- -v$E[flip, ]
    a <- limma::camera(v, idx, design, contrast = ci, inter.gene.cor = NA, sort = FALSE)
    b <- limma::camera(v, idx, design, contrast = ci, inter.gene.cor = 0.01, sort = FALSE)
    if (length(flip)) v$E[flip, ] <- -v$E[flip, ]
    out[[length(out) + 1]] <- data.table(set_key = k, n_genes = nrow(ss), direction = a$Direction, p = a$PValue,
                                         direction_preset = b$Direction, p_preset = b$PValue)
  }
  rbindlist(out)
}

fit_method <- function(counts, u, term, s = NULL) {
  d <- design_for(u, term)
  if (is.null(d)) return(NULL)
  design <- d$design
  y <- DGEList(counts); keep <- filterByExpr(y, design)
  if (sum(keep) < 50) return(NULL)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  coef <- grep(paste0("^", term), colnames(design), value = TRUE)[1]
  cam <- NULL
  if (method == "edger_ql") {
    y <- estimateDisp(y, design); fit <- glmQLFit(y, design, robust = TRUE)
    r <- topTags(glmQLFTest(fit, coef = coef), n = Inf, sort.by = "none")$table
    tt <- data.table(gene = rownames(r), logFC = r$logFC, p = r$PValue, padj = r$FDR, stat = r$F)
  } else if (method == "limma_voom") {
    v <- limma::voom(y, design)
    fit <- limma::eBayes(limma::lmFit(v, design), robust = TRUE)
    r <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
    tt <- data.table(gene = rownames(r), logFC = r$logFC, p = r$P.Value, padj = r$adj.P.Val, stat = r$t)
    if (!is.null(s) && nrow(s)) cam <- camera_sets(v, design, coef, s)
  } else if (method == "deseq2") {
    dds <- DESeqDataSetFromMatrix(round(y$counts), colData = data.frame(unit = colnames(y), row.names = colnames(y)),
                                  design = ~1)
    dds <- suppressMessages(DESeq(dds, full = design, betaPrior = FALSE, quiet = TRUE, parallel = ncpu > 1))
    r <- results(dds, contrast = as.numeric(colnames(design) == coef), independentFiltering = FALSE)
    tt <- data.table(gene = rownames(r), logFC = r$log2FoldChange, p = r$pvalue, padj = r$padj, stat = r$stat)
  } else {
    x <- cpm(y); lc <- cpm(y, log = TRUE)
    if (term == "group") {
      g <- d$u$group == "case"
      p <- apply(x, 1, function(z) suppressWarnings(wilcox.test(z[g], z[!g])$p.value))
      lfc <- rowMeans(lc[, g, drop = FALSE]) - rowMeans(lc[, !g, drop = FALSE])
    } else {
      n <- ncol(x); rho <- as.numeric(cor(apply(x, 1, rank), rank(d$u$stage)))  # units x genes ranks vs stage ranks
      p <- 2 * pt(-abs(rho * sqrt((n - 2) / pmax(1 - rho^2, 1e-12))), n - 2); lfc <- rho
    }
    tt <- data.table(gene = rownames(x), logFC = lfc, p = p, padj = p.adjust(p, "BH"), stat = NA_real_)
  }
  tt[, covariates := if (method == "wilcoxon") "" else paste(d$covs, collapse = "+")]
  list(tt = tt, cam = cam)
}

read_counts <- function(fl) {
  counts <- fread(cmd = paste("gzip -dc", shQuote(file.path(indir, fl))))
  gn <- counts[[1]]; counts <- as.matrix(counts[, -1, with = FALSE]); rownames(counts) <- gn
  counts
}
sets_for <- function(st, cn) if (is.null(sets)) NULL else sets[sets$stem == st & sets$contrast == cn]  # (names unlike the columns: data.table scoping)
safe <- function(expr, what) tryCatch(expr, error = function(e) { message(what, ": ", conditionMessage(e)); NULL })

# every fit: group contrasts with run == TRUE, and the stage per file when configured
fits <- unique(man[run == TRUE, .(region, cell_class, file, contrast)])
if (nzchar(stage_contrast)) fits <- rbind(fits, unique(man[, .(region, cell_class, file, contrast = stage_contrast)]))
fits <- unique(fits)

if (mode == "fit") {
  dir.create(file.path(outdir, paste0("de_", method)), recursive = TRUE, showWarnings = FALSE)
  log <- list(); cams <- list()
  for (fl in unique(fits$file)) {
    counts <- read_counts(fl); u <- units[match(colnames(counts), units$unit)]
    stem <- sub("\\.csv\\.gz$", "", fl)
    for (j in which(fits$file == fl)) {
      con <- fits$contrast[j]
      if (con == stage_contrast) {
        b <- stage_of(u); sel <- !is.na(b)
        if (sum(sel) < 6 || length(unique(b[sel])) < 3) next
        uu <- u[sel]; uu$stage <- b[sel]; term <- "stage"; nc <- sum(sel); nr <- NA
      } else {
        g <- u[[con]]; sel <- g %in% c("case", "ref")
        uu <- u[sel]; uu$group <- factor(g[sel], levels = c("ref", "case")); term <- "group"
        nc <- sum(uu$group == "case"); nr <- sum(uu$group == "ref")
      }
      r <- safe(fit_method(counts[, sel, drop = FALSE], uu, term, sets_for(stem, con)), paste(method, stem, con))
      if (is.null(r)) next
      fwrite(r$tt, file.path(outdir, paste0("de_", method), paste0(stem, "__", con, ".tsv.gz")), sep = "\t")
      log[[length(log) + 1]] <- data.table(file = fl, region = fits$region[j], cell_class = fits$cell_class[j],
        contrast = con, n_case = nc, n_ref = nr, n_genes = nrow(r$tt), n_fdr05 = sum(r$tt$padj < 0.05, na.rm = TRUE),
        covariates = r$tt$covariates[1])
      if (!is.null(r$cam) && nrow(r$cam)) cams[[length(cams) + 1]] <- cbind(data.table(stem = stem, contrast = con), r$cam)
    }
  }
  fwrite(rbindlist(log), file.path(outdir, paste0("de_log_", method, ".tsv")), sep = "\t")
  if (method == "limma_voom") fwrite(rbindlist(cams), file.path(outdir, "camera.tsv"), sep = "\t")
  message(basename(dsdir), " ", method, ": ", length(log), " fits, ", sum(vapply(cams, nrow, 1L)), " camera tests")
} else {
  set.seed(seed)
  dir.create(file.path(outdir, "null"), recursive = TRUE, showWarnings = FALSE)
  if (nrow(fits) > max_fits) fits <- fits[sort(sample.int(nrow(fits), max_fits))]
  null_sets <- if (is.null(sets)) NULL else unique(sets[, .(set_key, gene, sign)])
  keys <- if (is.null(null_sets)) character() else unique(null_sets$set_key)
  rows <- list(); hits <- list(); cams <- list()
  for (fl in unique(fits$file)) {
    counts <- read_counts(fl); u <- units[match(colnames(counts), units$unit)]
    stem <- sub("\\.csv\\.gz$", "", fl)
    for (j in which(fits$file == fl)) {
      con <- fits$contrast[j]
      if (con == stage_contrast) {
        b <- stage_of(u); sel <- !is.na(b)
        if (sum(sel) < 6 || length(unique(b[sel])) < 3) next
      } else {
        g <- u[[con]]; sel <- g %in% c("case", "ref")
        if (sum(g[sel] == "case") < 3 || sum(g[sel] == "ref") < 3) next
      }
      for (k in seq_len(perms)) {
        uu <- u[sel]
        if (con == stage_contrast) { uu$stage <- sample(b[sel]); term <- "stage"; nc <- sum(sel); nr <- NA
        } else { uu$group <- factor(sample(g[sel]), levels = c("ref", "case")); term <- "group"
                 nc <- sum(uu$group == "case"); nr <- sum(uu$group == "ref") }
        s <- if (length(keys)) null_sets[set_key %in% sample(keys, min(30, length(keys)))] else NULL
        r <- safe(fit_method(counts[, sel, drop = FALSE], uu, term, s), paste(method, stem, con, "perm", k))
        if (is.null(r)) next
        tt <- r$tt[!is.na(p)]
        chi <- qchisq(pmax(tt$p, 1e-300), 1, lower.tail = FALSE)
        rows[[length(rows) + 1]] <- data.table(file = fl, region = fits$region[j], cell_class = fits$cell_class[j],
          contrast = con, kind = term, perm = k, n_case = nc, n_ref = nr, n_genes = nrow(tt),
          frac_p05 = mean(tt$p < 0.05), n_fdr05 = sum(tt$padj < 0.05, na.rm = TRUE),
          lambda = median(chi) / qchisq(0.5, 1), covariates = tt$covariates[1])
        h <- tt[padj < 0.05, gene]
        if (length(h)) hits[[length(hits) + 1]] <- data.table(file = fl, contrast = con, perm = k, gene = h)
        if (!is.null(r$cam) && nrow(r$cam))
          cams[[length(cams) + 1]] <- cbind(data.table(stem = stem, contrast = con, kind = term, perm = k,
                                                       n_case = nc, n_ref = nr), r$cam)
      }
    }
  }
  res <- rbindlist(rows)
  fwrite(res, file.path(outdir, "null", paste0("de_null_", method, ".tsv")), sep = "\t")
  hh <- rbindlist(hits)
  if (!nrow(hh)) hh <- data.table(file = character(), contrast = character(), perm = integer(), gene = character())
  fwrite(hh, file.path(outdir, "null", paste0("fdr_hits_", method, ".tsv")), sep = "\t")
  if (method == "limma_voom") fwrite(rbindlist(cams), file.path(outdir, "null", "camera_null.tsv"), sep = "\t")
  if (nrow(res)) cat(basename(dsdir), method, nrow(res), "null fits; mean share of genes at p < 0.05:",
                     round(mean(res$frac_p05), 4), "\n")
}
