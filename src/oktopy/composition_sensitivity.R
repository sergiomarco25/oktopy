#!/usr/bin/env Rscript
# propeller (Phipson 2022: speckle's propeller.ttest with the logit transform), written with limma because speckle is
# not in the engine's R environment. Per unit (donor x region): each cell class's share with 0.5 added to every count,
# logit; lmFit ~ group (or ~ stage), eBayes(robust = TRUE) across classes; one test per class. No covariates, as in
# propeller and in the CLR primary. Sensitivity method of composition.case_control (phase 7b).
# Usage: Rscript composition_sensitivity.R IN_DIR OUT_DIR fit|null PERMS SEED [STAGE_CONTRAST]
#   IN_DIR/composition_counts.tsv  unit x class nuclei counts (oktopy.multiverse writes both inputs)
#   IN_DIR/composition_units.tsv   unit, region, one column per contrast (case / ref / empty), and the stage value in
#                                  the column named STAGE_CONTRAST
# fit   OUT_DIR/propeller.tsv
# null  labels (or stage values) shuffled among a region's units: OUT_DIR/null/composition_null_propeller.tsv
suppressPackageStartupMessages({ library(data.table); library(limma) })
args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]; outdir <- args[2]; mode <- args[3]; perms <- as.integer(args[4]); seed <- as.integer(args[5])
stage_contrast <- if (length(args) >= 6) args[6] else ""
stopifnot(mode %in% c("fit", "null"))
cnt <- fread(file.path(indir, "composition_counts.tsv"), colClasses = list(character = "unit"))
uts <- fread(file.path(indir, "composition_units.tsv"), colClasses = "character")
m <- as.matrix(cnt[, -1, with = FALSE]); rownames(m) <- cnt$unit
pp <- (m + 0.5) / rowSums(m + 0.5)
lt <- t(log(pp / (1 - pp)))  # classes x units
cons <- setdiff(names(uts), c("unit", "region"))

test <- function(units, x) {
  design <- model.matrix(~ x)
  fit <- eBayes(lmFit(lt[, units, drop = FALSE], design), robust = TRUE)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.table(cell_class = rownames(tt), effect = tt$logFC, t = tt$t, p = tt$P.Value)
}

set.seed(seed)
rows <- list()
for (reg in sort(unique(uts$region))) for (con in cons) {
  u <- uts[uts$region == reg]
  if (con == stage_contrast) {
    b <- suppressWarnings(as.numeric(u[[con]])); sel <- !is.na(b)
    if (sum(sel) < 6 || length(unique(b[sel])) < 3) next
    x0 <- b[sel]; nc <- sum(sel); nr <- 0L; mk <- function(x) x
  } else {
    g <- u[[con]]; sel <- g %in% c("case", "ref")
    if (sum(g[sel] == "case") < 3 || sum(g[sel] == "ref") < 3) next
    x0 <- g[sel]; nc <- sum(x0 == "case"); nr <- sum(x0 == "ref"); mk <- function(x) factor(x, levels = c("ref", "case"))
  }
  units <- u$unit[sel]
  for (k in if (mode == "fit") 0L else seq_len(perms)) {
    x <- if (k == 0L) x0 else sample(x0)
    r <- tryCatch(test(units, mk(x)), error = function(e) { message(reg, " ", con, ": ", conditionMessage(e)); NULL })
    if (!is.null(r)) rows[[length(rows) + 1]] <- cbind(data.table(region = reg, contrast = con, perm = k,
                                                                  n_case = nc, n_ref = nr), r)
  }
}
res <- rbindlist(rows)
if (!nrow(res)) res <- data.table(region = character(), contrast = character(), perm = integer(), n_case = integer(),
                                  n_ref = integer(), cell_class = character(), effect = numeric(), t = numeric(), p = numeric())
if (mode == "fit") {
  fwrite(res, file.path(outdir, "propeller.tsv"), sep = "\t")
} else {
  dir.create(file.path(outdir, "null"), recursive = TRUE, showWarnings = FALSE)
  fwrite(res, file.path(outdir, "null", "composition_null_propeller.tsv"), sep = "\t")
}
message("propeller ", mode, ": ", nrow(res), " tests", if (nrow(res) && mode == "null")
  paste0("; share p < 0.05 = ", round(mean(res$p < 0.05), 3)) else "")
