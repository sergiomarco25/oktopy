"""Method multiverse (phase 7b): claims re-tested with the sensitivity methods of the method registry
(oktopy/methods.toml) beside their primary, after each method has had the same shuffled-label check as the
primaries. Writes only to results/multiverse/; no verdict, store or atlas output changes (how disagreement between
methods counts is rating version 3, phase 7d).

  oktopy multiverse sets DATASET...
      the gene sets the donor score receives in each fit, recorded from a dry re-test run: <dataset>/sets.tsv (camera)
  oktopy multiverse fit METHOD DATASET...
      the method on the real labels: DE tables (de_<method>/), camera.tsv, or propeller.tsv
  oktopy multiverse null METHOD DATASET... [--perms N] [--max-fits 40] [--max-sets 30] [--seed 1]
      shuffled-label calibration of group contrasts and the stage: <dataset>/null/
  oktopy multiverse tests METHOD DATASET...
      internal and cross re-tests with METHOD in place of its task's primary: tests/<method>/<mode>/<dataset>.tsv
  oktopy multiverse all DATASET...
      every step above for each dataset, each in its own process (the workflow's multiverse rule runs this)
  oktopy multiverse summary
      calibration_summary.tsv, agreement.tsv (claim x dataset rows), agreement_summary.tsv, edger_refit_check.tsv

Methods (registry ids):
  edger_ql             the primary DE model refitted: fit must reproduce the DE cache; null adds the stage (edger_ql)
  limma_voom           voom + limma, single-gene claims (limma_voom); its fit and null also run camera
  deseq2               DESeq2 Wald test, single-gene claims (deseq2)
  wilcoxon             Wilcoxon rank-sum per gene, Spearman for the stage (wilcoxon_pseudobulk, spearman_pseudobulk)
  camera               competitive set test allowing for inter-gene correlation (camera); results come from limma_voom
  donor_score          the primary set test, null only: adds the stage to the phase 7a calibration (donor_score)
  donor_score_adjusted the donor score in a linear model with the DE cache's covariates (donor_score_adjusted)
  propeller            composition: logit shares, limma (propeller; Phipson 2022)
One dataset: `sets`, then `fit` and `null` of each method, then `tests` of each sensitivity method; then `summary`.
DE methods change single-gene claims only (sets of 2+ genes use the donor score, which does not depend on the DE
model); set methods change the donor-score claims; propeller changes the CLR composition claims.
"""
import argparse
import hashlib
import subprocess
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from oktopy import config as cfg
from oktopy import de_cache

HERE = Path(__file__).resolve().parent
DE = cfg.RESULTS / "de_cache"
OUT = cfg.RESULTS / "multiverse"
TESTS = cfg.RESULTS / "tests"
ALPHA = 0.05
BAND = (0.025, 0.075)
METHODS = {
    "edger_ql": "de", "limma_voom": "de", "deseq2": "de", "wilcoxon": "de",
    "camera": "set", "donor_score": "set", "donor_score_adjusted": "set", "propeller": "composition",
}
SENSITIVITY = ["limma_voom", "deseq2", "wilcoxon", "camera", "donor_score_adjusted", "propeller"]
LABELS = {"camera": "camera (limma-voom, inter-gene correlation estimated)",
          "donor_score_adjusted": "donor score, linear model with covariates"}


def set_key(genes_dir):
    """Identifier of a signed gene set (sign 0 is tested as up, as in retest.donor_score_test)."""
    items = sorted(f"{g}:{1 if s >= 0 else -1}" for g, s in genes_dir.items())
    return hashlib.sha1("|".join(items).encode()).hexdigest()[:16]


def applies(df, method):
    """Rows of a re-test table whose primary test the method replaces."""
    kind = METHODS[method]
    st = df.get("set_test", pd.Series("", index=df.index)).fillna("").astype(str)
    if kind == "de":
        return st.str.contains("single gene (FDR)", regex=False)
    if kind == "set":
        return st.str.contains("donor s", regex=False)
    return df.get("template_used", pd.Series("", index=df.index)).eq("composition_clr")


def method_ran(set_test, method):
    """Did the method itself produce this re-test row (not a fallback such as the random-gene set test)?"""
    st = "" if not isinstance(set_test, str) else set_test
    kind = METHODS[method]
    if kind == "de":
        return "single gene (FDR)" in st
    if kind == "composition":
        return "propeller" in st
    return LABELS[method] in st


def _stage_args():
    st = de_cache.stage_spec()
    code = de_cache.stage_contrast(st)
    return [code, ",".join(st["columns"]), str(st.get("group_column", "")), str(st.get("zero_group", ""))] if code else []


# ------------------------------------------------------------------------------------------------ inputs and fits

def collect_sets(ds):
    """Dry run of the internal and cross re-tests that records what the donor score is asked to test."""
    from oktopy import retest as rt
    rec = []

    def recorder(dataset, region, cls, contrast, genes_dir, ut=None):
        rec.append((region, cls, contrast, dict(genes_dir)))
        return None

    rt.donor_score_test = recorder
    rt.state_existence_test = lambda *a, **k: None
    for mode in ("internal", "cross"):
        rt.run(mode, ds, out_dir=OUT / ds / "_dry")
    rows = [{"stem": f"{reg}__{rt.fname(cls)}", "contrast": con, "set_key": set_key(gd), "gene": g,
             "sign": 1 if s >= 0 else -1} for reg, cls, con, gd in rec for g, s in gd.items()]
    df = pd.DataFrame(rows, columns=["stem", "contrast", "set_key", "gene", "sign"]).drop_duplicates()
    (OUT / ds).mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / ds / "sets.tsv", sep="\t", index=False)
    print(ds, "sets:", df.set_key.nunique(), "gene sets in", df[["stem", "contrast"]].drop_duplicates().shape[0], "fits")
    return 0


def run_de(ds, method, mode, perms=0, max_fits=0, seed=0):
    sets = OUT / ds / "sets.tsv"
    sets = str(sets) if method == "limma_voom" and sets.is_file() else "-"
    cmd = [de_cache.rscript(), "--vanilla", str(HERE / "de_sensitivity.R"), method, str(DE / ds), str(OUT / ds), sets,
           mode, str(perms), str(max_fits), str(seed)] + _stage_args()
    (OUT / ds).mkdir(parents=True, exist_ok=True)
    return subprocess.run(cmd, check=False).returncode


def write_composition_inputs(ds):
    from oktopy import retest as rt
    from oktopy.contrasts import CONTRASTS, stage_value
    o = rt.cell_obs(ds)
    o = o.assign(unit=o.donor.astype(str) + "|" + o.region.astype(str))
    counts = pd.crosstab(o.unit, o.cell_class.astype(str))
    meta = rt.unit_meta(ds)
    units = meta[["region"] + [c for c in CONTRASTS if c in meta]].copy()
    units["region"] = units.region.astype(str)
    code = de_cache.stage_contrast()
    if code:
        units[code] = stage_value(meta)
    keep = counts.index.intersection(units.index)
    d = OUT / ds / "composition"
    d.mkdir(parents=True, exist_ok=True)
    counts.loc[keep].rename_axis("unit").to_csv(d / "composition_counts.tsv", sep="\t")
    units.loc[keep].rename_axis("unit").to_csv(d / "composition_units.tsv", sep="\t")
    return d


def run_propeller(ds, mode, perms=0, seed=0):
    d = write_composition_inputs(ds)
    cmd = [de_cache.rscript(), "--vanilla", str(HERE / "composition_sensitivity.R"), str(d), str(OUT / ds), mode,
           str(perms), str(seed), de_cache.stage_contrast()]
    return subprocess.run(cmd, check=False).returncode


# ------------------------------------------------------------------------------------------------ set tests

def _covariates(u):
    """The DE cache's covariate rule (de_cache.R fit_one): age, sex, PMI when known for every unit and not constant."""
    out = []
    age = pd.to_numeric(u.get("age", pd.Series(np.nan, index=u.index)), errors="coerce")
    if age.notna().all() and age.std() > 0:
        out.append(("age_s", ((age - age.mean()) / age.std()).values))
    sex = u.get("sex", pd.Series("", index=u.index)).astype(str)
    if sex.isin(["M", "F"]).all() and sex.nunique() > 1:
        out.append(("sex", sex.eq("M").astype(float).values))
    pmi = pd.to_numeric(u.get("pmi", pd.Series(np.nan, index=u.index)), errors="coerce")
    if pmi.notna().all() and pmi.std() > 0:
        out.append(("pmi_s", ((pmi - pmi.mean()) / pmi.std()).values))
    return out


def ols_fit(y, x, u):
    """y ~ x + covariates (dropped from the end until the design is estimable with >= 2 residual df, as in
    de_cache.R): x's coefficient, its standard error and two-sided p, the residual SD and the covariates kept.
    None if the model cannot be fitted."""
    covs = _covariates(u)
    n = len(y)
    while True:
        X = np.column_stack([np.ones(n), x] + [c for _, c in covs])
        if np.linalg.matrix_rank(X) == X.shape[1] and n - X.shape[1] >= 2:
            break
        if not covs:
            return None
        covs = covs[:-1]
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    df = n - X.shape[1]
    resid = y - X @ beta
    sigma = np.sqrt(resid @ resid / df)
    se = sigma * np.sqrt(np.linalg.inv(X.T @ X)[1, 1])
    if not np.isfinite(se) or se == 0:
        return None
    return {"beta": float(beta[1]), "se": float(se), "p": float(2 * stats.t.sf(abs(beta[1] / se), df)),
            "sigma": float(sigma), "covs": [c for c, _ in covs]}


def ols_term(y, x, u):
    """Coefficient and two-sided p of x in y ~ x + covariates, and the covariates kept (ols_fit)."""
    f = ols_fit(y, x, u)
    return (None, None, []) if f is None else (f["beta"], f["p"], f["covs"])


def native_smd(fit, x, stage):
    """The adjusted coefficient as a standardised effect on the scale rating.to_smd gives unadjusted tests, with its
    own standard error, so a covariate that is collinear with x widens the interval instead of shrinking the effect.
    Two groups: d = beta / residual SD, se as in to_smd. Stage: d = 2 b / sqrt(1 - r^2) from the standardised slope
    b = beta sd(x) / residual SD (r = b / sqrt(1 + b^2), the correlation of an unadjusted fit), r capped at 0.95."""
    if fit["sigma"] <= 0:
        return np.nan, np.nan
    if not stage:
        d = float(np.clip(fit["beta"] / fit["sigma"], -3, 3))
        n = len(x)
        return d, float(np.sqrt((fit["se"] / fit["sigma"]) ** 2 + d ** 2 / (2 * n)))
    sx = float(np.std(x, ddof=1))
    b = fit["beta"] * sx / fit["sigma"]
    b = float(np.clip(b, -0.95 / np.sqrt(1 - 0.95 ** 2), 0.95 / np.sqrt(1 - 0.95 ** 2)))
    return 2 * b, float(2 * fit["se"] * sx / fit["sigma"] * np.sqrt(1 + b ** 2))


def adjusted_score_test(dataset, region, cls, contrast, genes_dir, ut=None):
    """retest.donor_score_test's score in a linear model with the DE cache's covariates (t-test of the group or stage
    coefficient) instead of Mann-Whitney / Spearman."""
    from oktopy import retest as rt
    from oktopy.contrasts import stage_value
    zt, ut = rt.unit_logcpm_z(dataset, region, cls), (rt.unit_table(dataset) if ut is None else ut)
    if zt is None or ut is None:
        return None
    hit, up = {}, rt.gene_index(zt.index)
    for g, sgn in genes_dir.items():
        m = rt.map_gene(g, up)
        if m is not None:
            hit[m] = sgn if sgn != 0 else 1
    if len(hit) < rt.MIN_SCORE_GENES:
        return None
    score = zt.loc[list(hit)].mul(pd.Series(hit), axis=0).mean(0)
    if contrast == rt.STAGE_CONTRAST:
        b = stage_value(ut.loc[score.index]).dropna()
        if len(b) < 6 or b.nunique() < 3:
            return None
        units, x, n_case, n_ref = b.index, b.astype(float).values, len(b), 0
        label = f"{LABELS['donor_score_adjusted']} ~ {rt.STAGE_LABEL}"
    else:
        if contrast not in ut:
            return None
        g = ut.loc[score.index, contrast].fillna("")
        units = g.index[g.isin(["case", "ref"])]
        x = g.loc[units].eq("case").astype(float).values
        n_case, n_ref = int(x.sum()), int(len(x) - x.sum())
        if n_case < 3 or n_ref < 3:
            return None
        label = LABELS["donor_score_adjusted"]
    fit = ols_fit(score.loc[units].values, x, ut.loc[units])
    if fit is None:
        return None
    beta, p = fit["beta"], fit["p"]
    d, se_d = native_smd(fit, x, contrast == rt.STAGE_CONTRAST)
    return {"z": np.sign(beta) * stats.norm.isf(max(p, 1e-300) / 2), "p": p, "effect": beta, "n_case": n_case,
            "n_ref": n_ref, "n_genes": len(hit), "test": label, "d_native": d, "se_native": se_d}


@lru_cache(maxsize=256)
def method_table(dataset, method, stem_con):
    """A sensitivity method's DE table in the DE cache's format (gene index; logFC, p, padj, z)."""
    f = OUT / dataset / f"de_{method}" / f"{stem_con}.tsv.gz"
    if not f.exists():
        return None
    t = pd.read_csv(f, sep="\t").dropna(subset=["p"]).set_index("gene")
    t["padj"] = t.padj.fillna(1.0)
    t["z"] = np.sign(t.logFC) * stats.norm.isf(np.clip(t.p, 1e-300, 1) / 2)
    return t


def camera_test_for(ds):
    from oktopy import retest as rt
    f = OUT / ds / "camera.tsv"
    cam = pd.read_csv(f, sep="\t") if f.is_file() and f.stat().st_size > 5 else pd.DataFrame(
        columns=["stem", "contrast", "set_key", "n_genes", "direction", "p"])
    cam = cam.drop_duplicates(["stem", "contrast", "set_key"]).set_index(["stem", "contrast", "set_key"])
    lf = OUT / ds / "de_log_limma_voom.tsv"
    log = pd.read_csv(lf, sep="\t") if lf.is_file() and lf.stat().st_size > 5 else pd.DataFrame(
        columns=["region", "cell_class", "contrast", "n_case", "n_ref"])

    def test(dataset, region, cls, contrast, genes_dir, ut=None):
        stem = f"{region}__{rt.fname(cls)}"
        key = (stem, contrast, set_key(genes_dir))
        if ut is not None or key not in cam.index:
            return None
        r = cam.loc[key]
        p = float(r.p)
        t = method_table(dataset, "limma_voom", f"{stem}__{contrast}")
        g = [x for x in genes_dir if t is not None and x in t.index]
        eff = float(np.mean([t.loc[x, "logFC"] * (1 if genes_dir[x] >= 0 else -1) for x in g])) if g else np.nan
        n = log[(log.region == region) & (log.cell_class == cls) & (log.contrast == contrast)]
        return {"z": (1 if r.direction == "Up" else -1) * stats.norm.isf(max(p, 1e-300) / 2), "p": p, "effect": eff,
                "n_case": int(n.n_case.iloc[0]) if len(n) else np.nan,
                "n_ref": int(n.n_ref.iloc[0]) if len(n) and pd.notna(n.n_ref.iloc[0]) else 0,
                "n_genes": int(r.n_genes), "test": LABELS["camera"]}
    return test


def propeller_test_for(ds, orig):
    from oktopy import retest as rt
    f = OUT / ds / "propeller.tsv"
    tab = pd.read_csv(f, sep="\t", dtype={"region": str}) if f.is_file() and f.stat().st_size > 5 else None

    def test(dataset, regions, classes, contrast, expected, units=None, meta=None):
        if units is not None or meta is not None:
            return orig(dataset, regions, classes, contrast, expected, units=units, meta=meta)
        if tab is None:
            return pd.DataFrame()
        rows = []
        for reg in regions:
            for cls in classes:
                r = tab[(tab.region == reg) & (tab.cell_class == cls) & (tab.contrast == contrast)]
                if not len(r):
                    continue
                r = r.iloc[0]
                rows.append({"region": reg, "cell_class": cls, "p": float(r.p), "effect": float(r.effect),
                             "z": np.sign(r.effect) * stats.norm.isf(max(float(r.p), 1e-300) / 2),
                             "n_case": int(r.n_case), "n_ref": int(r.n_ref),
                             "test": "propeller (logit shares, limma)"})
        return pd.DataFrame(rows)
    return test


def run_tests(ds, method):
    """Internal and cross re-tests with the method in place of its task's primary. State-existence and genetic tests
    are skipped (no method here changes them; the field's genetic hook may read edgeR-only columns)."""
    from oktopy import retest as rt
    rt.state_existence_test = lambda *a, **k: None
    rt.genetic_test = lambda *a, **k: []
    kind = METHODS[method]
    if kind == "de":
        rt.de_table = lambda dataset, region, cls, contrast: method_table(
            dataset, method, f"{region}__{rt.fname(cls)}__{contrast}")
    elif method == "camera":
        rt.donor_score_test = camera_test_for(ds)
    elif method == "donor_score_adjusted":
        rt.donor_score_test = adjusted_score_test
    elif method == "propeller":
        rt.composition_test = propeller_test_for(ds, rt.composition_test)
    else:
        print(method, "is not a sensitivity method (no tests)")
        return 2
    for mode in ("internal", "cross"):
        rt.run(mode, ds, out_dir=OUT / "tests" / method)
    return 0


# ------------------------------------------------------------------------------------------------ calibration

def null_scores(ds, method, perms, max_fits, max_sets, seed):
    """Shuffled-label check of the donor score (primary) or the adjusted donor score, group contrasts and the stage."""
    from oktopy import calibration as cal
    from oktopy import retest as rt
    from oktopy.contrasts import stage_value
    test = rt.donor_score_test if method == "donor_score" else adjusted_score_test
    rng = np.random.default_rng(seed)
    ut, log = rt.unit_table(ds), rt.de_log(ds)
    if ut is None or log.empty:
        print(ds, "no unit table or DE fits: nothing to calibrate")
        return 0
    combos = log[["region", "cell_class", "contrast"]].drop_duplicates()
    combos = combos[combos.contrast.isin(ut.columns) | combos.contrast.eq(rt.STAGE_CONTRAST)].reset_index(drop=True)
    if len(combos) > max_fits:
        combos = combos.loc[sorted(rng.choice(len(combos), max_fits, replace=False))]
    sets = cal.claim_gene_sets(rt, 2)
    if len(sets) > max_sets:
        sets = [sets[i] for i in sorted(rng.choice(len(sets), max_sets, replace=False))]
    stage_col = next((c for c in de_cache.stage_spec().get("columns", []) if c in ut), None)
    rows = []
    for reg, cls, con in combos.itertuples(index=False):
        zt = rt.unit_logcpm_z(ds, reg, cls)
        if zt is None:
            continue
        idx = [u for u in zt.columns if u in ut.index]
        stage = con == rt.STAGE_CONTRAST
        if stage:
            if stage_col is None:
                continue
            b = stage_value(ut.loc[idx]).dropna()
            if len(b) < 6 or b.nunique() < 3:
                continue
        for k in range(perms):
            up = ut.copy()
            if stage:
                up[stage_col] = up[stage_col].astype(object)
                up.loc[idx, stage_col] = np.nan
                up.loc[b.index, stage_col] = rng.permutation(b.values)
            else:
                up.loc[idx, con] = rng.permutation(ut.loc[idx, con].fillna("").values)
            for label, signs in sets:
                r = test(ds, reg, cls, con, signs, ut=up)
                if r is not None:
                    rows.append({"region": reg, "cell_class": cls, "contrast": con, "kind": "stage" if stage else "group",
                                 "perm": k + 1, "set": label, "n_genes": r["n_genes"], "n_case": r["n_case"],
                                 "n_ref": r["n_ref"], "p": r["p"]})
    res = pd.DataFrame(rows)
    (OUT / ds / "null").mkdir(parents=True, exist_ok=True)
    res.to_csv(OUT / ds / "null" / f"scores_null_{method}.tsv", sep="\t", index=False)
    if len(res):
        print(ds, method, len(res), "null tests; share p < 0.05:",
              res.groupby("kind").p.apply(lambda p: round(float((p < ALPHA).mean()), 3)).to_dict())
    return 0


def _status(share):
    return "calibrated" if BAND[0] <= share <= BAND[1] else ("liberal" if share > BAND[1] else "conservative")


def _read(f):
    return pd.read_csv(f, sep="\t") if f.is_file() and f.stat().st_size > 5 else None


def calibration_rows(ds):
    from oktopy.calibration import wilson
    d, rows = OUT / ds / "null", []

    def add(method, test, kind, p, units, extra=None):
        k, n = int((p < ALPHA).sum()), len(p)
        if not n:
            return
        lo, hi = wilson(k, n)
        rows.append({"dataset": ds, "method": method, "test": test, "kind": kind, "n_null": n, "share_p05": k / n,
                     "ci_lo": lo, "ci_hi": hi, "median_n_units": float(np.nanmedian(units)) if len(units) else np.nan,
                     "status": _status(k / n), **(extra or {})})

    for m in ("edger_ql", "limma_voom", "deseq2", "wilcoxon"):
        t = _read(d / f"de_null_{m}.tsv")
        if t is None:
            continue
        hits = _read(d / f"fdr_hits_{m}.tsv")
        for kind, part in t.groupby("kind"):
            share = float(part.frac_p05.mean())
            h = hits[hits.contrast.isin(part.contrast)] if hits is not None else None
            per_gene = float(h.groupby("gene").size().max() / len(part)) if h is not None and len(h) else 0.0
            rows.append({"dataset": ds, "method": m, "test": "genes", "kind": kind, "n_null": len(part),
                         "share_p05": share, "ci_lo": np.nan, "ci_hi": np.nan,
                         "median_n_units": float((part.n_case + part.n_ref.fillna(0)).median()),
                         "status": _status(share), "median_lambda": float(part["lambda"].median()),
                         "fits_with_fdr05_hit": float((part.n_fdr05 > 0).mean()), "max_gene_fdr05_share": per_gene})
    cam = _read(d / "camera_null.tsv")
    if cam is not None:
        for kind, part in cam.groupby("kind"):
            for size, sub in (("2-4 genes", part[part.n_genes <= 4]), ("5+ genes", part[part.n_genes >= 5])):
                units = sub.n_case + sub.n_ref.fillna(0)
                add("camera", size, kind, sub.p, units)
                add("camera_preset_cor", size, kind, sub.p_preset, units)
    for m in ("donor_score", "donor_score_adjusted"):
        t = _read(d / f"scores_null_{m}.tsv")
        if t is None:
            continue
        for kind, part in t.groupby("kind"):
            for size, sub in (("2-4 genes", part[part.n_genes <= 4]), ("5+ genes", part[part.n_genes >= 5])):
                add(m, size, kind, sub.p, sub.n_case + sub.n_ref)
    t = _read(d / "composition_null_propeller.tsv")
    if t is not None:
        for kind, part in t.assign(kind=np.where(t.n_ref.eq(0), "stage", "group")).groupby("kind"):
            add("propeller", "classes", kind, part.p, part.n_case + part.n_ref)
    return rows


# ------------------------------------------------------------------------------------------------ agreement

def agreement():
    rows = []
    for method in SENSITIVITY:
        for f in sorted((OUT / "tests" / method).glob("*/*.tsv")):
            prim_f = TESTS / f.parent.name / f.name
            if not prim_f.is_file():
                continue
            p, m = pd.read_csv(prim_f, sep="\t"), pd.read_csv(f, sep="\t")
            if p.empty or m.empty:  # e.g. no claim has its source data in this dataset (header-only internal file)
                continue
            key =[c for c in ["claim_id", "dataset", "contrast", "variant", "template_used"] if c in p and c in m]
            for t in (p, m):
                t[key] = t[key].fillna("")
                t["_k"] = t.groupby(key).cumcount()
            p = p[applies(p, method) & p.outcome.ne("not_testable")]
            j = p.merge(m[key + ["_k", "outcome", "p", "z", "set_test"]], on=key + ["_k"], how="left",
                        suffixes=("_primary", "_method"))
            for r in j.itertuples(index=False):
                po = r.outcome_primary
                mo = r.outcome_method if method_ran(r.set_test_method, method) else "not_testable"
                rows.append({"method": method, "mode": f.parent.name, "claim_id": r.claim_id, "dataset": r.dataset,
                             "contrast": r.contrast, "variant": getattr(r, "variant", ""),
                             "template_used": r.template_used, "primary_test": r.set_test_primary,
                             "method_test": r.set_test_method, "outcome_primary": po, "outcome_method": mo,
                             "p_primary": r.p_primary, "p_method": r.p_method, "z_primary": r.z_primary,
                             "z_method": r.z_method, "n_genes": getattr(r, "n_genes", np.nan),
                             "same_outcome": po == mo})
    return pd.DataFrame(rows)


def edger_check():
    """The edgeR refit (de_sensitivity.R, edger_ql) against the DE cache: largest differences per dataset."""
    rows = []
    for d in sorted(p for p in OUT.iterdir() if (p / "de_edger_ql").is_dir()) if OUT.is_dir() else []:
        n, dp, dl, missing = 0, 0.0, 0.0, 0
        for f in sorted((d / "de_edger_ql").glob("*.tsv.gz")):
            c = DE / d.name / "de" / f.name
            if not c.is_file():
                missing += 1
                continue
            a = pd.read_csv(c, sep="\t").set_index("gene")
            b = pd.read_csv(f, sep="\t").set_index("gene").reindex(a.index)
            n += 1
            dp = max(dp, float(np.nanmax(np.abs(a.p - b.p))) if b.p.notna().all() else np.inf)
            dl = max(dl, float(np.nanmax(np.abs(a.logFC - b.logFC))) if b.logFC.notna().all() else np.inf)
        rows.append({"dataset": d.name, "fits_compared": n, "refits_without_cache_table": missing,
                     "max_abs_diff_p": dp, "max_abs_diff_logFC": dl})
    return pd.DataFrame(rows)


def summary():
    cal = pd.DataFrame([r for d in sorted(p for p in OUT.iterdir() if p.is_dir() and p.name != "tests")
                        for r in calibration_rows(d.name)]) if OUT.is_dir() else pd.DataFrame()
    OUT.mkdir(parents=True, exist_ok=True)
    cal.round(4).to_csv(OUT / "calibration_summary.tsv", sep="\t", index=False)
    ag = agreement()
    if len(ag) and len(cal):  # did the method pass the shuffled-label check in this dataset, for this kind of test?
        ok = cal[cal.status.eq("calibrated")]
        kind = np.where(ag.contrast.astype(str).eq(de_cache.stage_contrast()), "stage", "group")
        test = np.select([ag.method.map(METHODS).eq("de"), ag.method.eq("propeller"), ag.n_genes.fillna(0) <= 4.5],
                         ["genes", "classes", "2-4 genes"], "5+ genes")
        keys = set(zip(ok.dataset, ok.method, ok.test, ok.kind))
        ag["method_calibrated"] = [(d, m, t, k) in keys for d, m, t, k in zip(ag.dataset, ag.method, test, kind)]
    ag.to_csv(OUT / "agreement.tsv", sep="\t", index=False)
    if len(ag):
        testable = ag[ag.outcome_method.ne("not_testable")]
        s = testable.groupby(["method", "dataset"]).agg(
            compared=("same_outcome", "size"), same_outcome=("same_outcome", "mean"),
            primary_reproduced=("outcome_primary", lambda o: int(o.eq("reproduced").sum()))).reset_index()
        lost = testable[testable.outcome_primary.eq("reproduced") & testable.outcome_method.ne("reproduced")]
        gained = testable[testable.outcome_primary.ne("reproduced") & testable.outcome_method.eq("reproduced")]
        s["lost"] = s.set_index(["method", "dataset"]).index.map(lost.groupby(["method", "dataset"]).size()).fillna(0).astype(int)
        s["gained"] = s.set_index(["method", "dataset"]).index.map(gained.groupby(["method", "dataset"]).size()).fillna(0).astype(int)
        s["not_testable_by_method"] = s.set_index(["method", "dataset"]).index.map(
            ag[ag.outcome_method.eq("not_testable")].groupby(["method", "dataset"]).size()).fillna(0).astype(int)
        s.round(3).to_csv(OUT / "agreement_summary.tsv", sep="\t", index=False)
        print(s.round(3).to_string(index=False))
    if len(cal):
        print(cal.drop(columns=[c for c in ("ci_lo", "ci_hi") if c in cal]).round(3).to_string(index=False))
    chk = edger_check()
    chk.to_csv(OUT / "edger_refit_check.tsv", sep="\t", index=False)
    if len(chk):
        print(chk.to_string(index=False))
    print("wrote", OUT / "calibration_summary.tsv", OUT / "agreement.tsv", OUT / "agreement_summary.tsv",
          OUT / "edger_refit_check.tsv")
    return 0


def run_all(ds, perms_de=5, perms_scores=10, perms_comp=100, max_fits=40, max_sets=30, seed=1):
    """Every step for one dataset, each in its own process (tests swap one engine function per method)."""
    cli = [sys.executable, "-m", "oktopy", "--project", str(cfg.REPO), "multiverse"]
    fits = ["--max-fits", str(max_fits), "--seed", str(seed)]
    steps = [["sets", ds]]
    for m in ("edger_ql", "limma_voom", "deseq2", "wilcoxon"):
        steps += [["fit", m, ds], ["null", m, ds, "--perms", str(perms_de), *fits]]
    steps += [["fit", "propeller", ds], ["null", "propeller", ds, "--perms", str(perms_comp), "--seed", str(seed)]]
    for m in ("donor_score", "donor_score_adjusted"):
        steps += [["null", m, ds, "--perms", str(perms_scores), "--max-sets", str(max_sets), *fits]]
    steps += [["tests", m, ds] for m in SENSITIVITY]
    for st in steps:
        print(">>", " ".join(st), flush=True)
        rc = subprocess.run(cli + st, check=False).returncode
        if rc:
            print("failed:", " ".join(st), "exit", rc, flush=True)
            return rc
    (OUT / ds / "done").touch()
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="oktopy multiverse", description=__doc__.split("\n\n")[0])
    ap.add_argument("what", choices=["sets", "fit", "null", "tests", "all", "summary"])
    ap.add_argument("args", nargs="*", help="METHOD DATASET... (fit, null, tests) or DATASET... (sets, all)")
    ap.add_argument("--perms", type=int)
    ap.add_argument("--max-fits", type=int, default=40)
    ap.add_argument("--max-sets", type=int, default=30)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(argv)
    if a.what == "summary":
        return summary()
    if a.what in ("sets", "all") and not a.args:
        ap.error("name at least one dataset")
    if a.what == "sets":
        return max(collect_sets(ds) for ds in a.args)
    if a.what == "all":
        return max(run_all(ds, max_fits=a.max_fits, max_sets=a.max_sets, seed=a.seed) for ds in a.args)
    if len(a.args) < 2 or a.args[0] not in METHODS:
        ap.error(f"give METHOD DATASET...; methods: {', '.join(METHODS)}")
    method, dss = a.args[0], a.args[1:]
    codes = []
    for ds in dss:
        if a.what == "tests":
            codes.append(run_tests(ds, method))
        elif method == "propeller":
            codes.append(run_propeller(ds, a.what, a.perms or 100, a.seed))
        elif METHODS[method] == "de" or method == "camera":
            m = "limma_voom" if method == "camera" else method
            codes.append(run_de(ds, m, a.what, a.perms or 5, a.max_fits, a.seed))
        elif a.what == "null":
            codes.append(null_scores(ds, method, a.perms or 10, a.max_fits, a.max_sets, a.seed))
        else:
            print(method, "needs no fit (computed in `tests`)")
            codes.append(0)
    return max(codes)


if __name__ == "__main__":
    sys.exit(main())
