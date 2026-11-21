"""Shuffled-label calibration of the engine's current tests (phase 7a; no verdict or store output changes).

Case and reference labels are shuffled among a fit's units, so no difference is real: a calibrated test then gives
p < 0.05 in about 5% of tests. Results go to results/robustness/calibration/<dataset>/.

  oktopy calibrate de DATASET [--perms 5] [--max-fits 40] [--seed 1]
      edgeR (de_calibration.R, the DE cache's model) on shuffled labels: de_null.tsv, null/*.tsv.gz
  oktopy calibrate sets DATASET [--max-sets 60] [--seed 1]
      the random-gene competitive test (claims' gene sets of 2-4 genes, and of 5+ for comparison) on the null DE
      tables: sets_null.tsv
  oktopy calibrate scores DATASET [--small] [--perms 20] [--max-fits 40] [--max-sets 30]
      the donor score test with shuffled unit labels on claims' gene sets of 5+ genes (scores_null.tsv), or with
      --small of 2-4 genes (scores_small_null.tsv)
  oktopy calibrate composition DATASET [--perms 100] [--seed 1]
      the CLR composition test with shuffled unit labels: composition_null.tsv
  oktopy calibrate state DATASET [--max-sets 40] [--seed 1]
      the state-existence test on random genes in the cell classes the state claims name: state_null.tsv (no labels
      to shuffle: random genes are its null, so p < 0.05 in about 5% and all three rules passed in at most 5%)
  oktopy calibrate summary
      one row per dataset x test: results/robustness/calibration_summary.tsv

Status: calibrated when the share of p < 0.05 lies in [0.025, 0.075], liberal above, conservative below. Tests on
one shuffle share units, so the Wilson interval reported beside the share is optimistic.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from oktopy import config as cfg
from oktopy import de_cache

HERE = Path(__file__).resolve().parent
DE = cfg.RESULTS / "de_cache"
OUT = cfg.RESULTS / "robustness" / "calibration"
BAND = (0.025, 0.075)
ALPHA = 0.05


def _null_table(f):
    t = pd.read_csv(f, sep="\t").dropna(subset=["p"]).set_index("gene")
    t["z"] = np.sign(t.logFC) * stats.norm.isf(np.clip(t.p, 1e-300, 1) / 2)
    o = np.argsort(t.p.values)
    q = t.p.values[o] * len(t) / np.arange(1, len(t) + 1)
    padj = np.empty(len(t))
    padj[o] = np.minimum.accumulate(q[::-1])[::-1].clip(max=1)
    t["padj"] = padj
    return t


def claim_gene_sets(rt, min_genes, max_genes=None):
    """(label, gene -> expected sign) for every claim signature variant, oriented as claimed (+1 for 'up')."""
    sig = rt.signatures()
    if sig is None or sig.empty:
        return []
    out = []
    keys = ["signature_id", "variant"] if "variant" in sig else ["signature_id"]
    for key, sv in sig.groupby(keys):
        signs = rt.expected_signs(sv, 1)
        if len(signs) >= min_genes and (max_genes is None or len(signs) <= max_genes):
            out.append(("|".join(map(str, key if isinstance(key, tuple) else (key,))), signs))
    return out


def run_de(ds, perms, max_fits, seed):
    cmd = [de_cache.rscript(), "--vanilla", str(HERE / "de_calibration.R"), str(DE / ds), str(OUT / ds),
           str(perms), str(max_fits), str(seed)]
    return subprocess.run(cmd, check=False).returncode


def run_sets(ds, max_sets, seed):
    from oktopy import retest as rt
    rng = np.random.default_rng(seed)
    files = sorted((OUT / ds / "null").glob("*.tsv.gz"))
    if not files:
        print(ds, "no null DE tables (no fits to calibrate, or `calibrate de` not run)")
        return 0
    sets = claim_gene_sets(rt, 2)
    if len(sets) > max_sets:
        sets = [sets[i] for i in sorted(rng.choice(len(sets), max_sets, replace=False))]
    rows = []
    for f in files:
        t = _null_table(f)
        stem = f.name[:-len(".tsv.gz")]
        up = rt.gene_index(t.index)
        for label, signs in sets:
            mapped = {}
            for g, s in signs.items():
                m = rt.map_gene(g, up)
                if m is not None:
                    mapped[m] = s
            if len(mapped) < 2:
                continue
            r = rt.gene_set_test(t, mapped)
            if r is not None:
                rows.append({"null_fit": stem, "set": label, "n_genes": r["n_genes"], "p": r["p"]})
    res = pd.DataFrame(rows)
    res.to_csv(OUT / ds / "sets_null.tsv", sep="\t", index=False)
    if len(res):
        for name, part in (("2-4 genes", res[res.n_genes <= 4]), ("5+ genes", res[res.n_genes >= 5])):
            if len(part):
                print(ds, f"competitive set test, {name}: {len(part)} null tests, share p < 0.05 = "
                      f"{(part.p < ALPHA).mean():.3f}")
    return 0


def run_scores(ds, perms, max_fits, max_sets, seed, small=False):
    from oktopy import retest as rt
    rng = np.random.default_rng(seed)
    ut = rt.unit_table(ds)
    log = rt.de_log(ds)
    if ut is None or log.empty:
        print(ds, "no unit table or DE fits: nothing to calibrate")
        return 0
    log = log[log.contrast.isin(ut.columns) & log.n_ref.notna()]
    combos = log[["region", "cell_class", "contrast"]].drop_duplicates().reset_index(drop=True)
    if len(combos) > max_fits:
        combos = combos.loc[sorted(rng.choice(len(combos), max_fits, replace=False))]
    sets = claim_gene_sets(rt, 2, 4) if small else claim_gene_sets(rt, 5)
    if len(sets) > max_sets:
        sets = [sets[i] for i in sorted(rng.choice(len(sets), max_sets, replace=False))]
    rows = []
    for reg, cls, con in combos.itertuples(index=False):
        zt = rt.unit_logcpm_z(ds, reg, cls)
        if zt is None:
            continue
        idx = [u for u in zt.columns if u in ut.index]
        for k in range(perms):
            up = ut.copy()
            up.loc[idx, con] = rng.permutation(ut.loc[idx, con].fillna("").values)
            for label, signs in sets:
                r = rt.donor_score_test(ds, reg, cls, con, signs, ut=up)
                if r is not None:
                    rows.append({"region": reg, "cell_class": cls, "contrast": con, "perm": k + 1, "set": label,
                                 "n_genes": r["n_genes"], "n_case": r["n_case"], "n_ref": r["n_ref"], "p": r["p"]})
    res = pd.DataFrame(rows)
    (OUT / ds).mkdir(parents=True, exist_ok=True)
    res.to_csv(OUT / ds / ("scores_small_null.tsv" if small else "scores_null.tsv"), sep="\t", index=False)
    if len(res):
        print(ds, f"donor score test{', 2-4 genes' if small else ''}: {len(res)} null tests, share p < 0.05 = "
              f"{(res.p < ALPHA).mean():.3f}")
    return 0


def run_composition(ds, perms, seed):
    from oktopy import retest as rt
    from oktopy.contrasts import CONTRASTS
    rng = np.random.default_rng(seed)
    meta = rt.unit_meta(ds)
    classes = sorted(rt.cell_obs(ds).cell_class.astype(str).unique())
    rows = []
    for con in [c for c in CONTRASTS if c in meta]:
        for reg in sorted(meta.region.astype(str).unique()):
            idx = meta.index[meta.region.astype(str).eq(reg) & meta[con].isin(["case", "ref"])]
            if (meta.loc[idx, con] == "case").sum() < 3 or (meta.loc[idx, con] == "ref").sum() < 3:
                continue
            for k in range(perms):
                mp = meta.copy()
                mp.loc[idx, con] = rng.permutation(meta.loc[idx, con].values)
                df = rt.composition_test(ds, [reg], classes, con, 0, meta=mp)
                for r in df.itertuples(index=False):
                    rows.append({"region": reg, "contrast": con, "perm": k + 1, "cell_class": r.cell_class,
                                 "n_case": r.n_case, "n_ref": r.n_ref, "p": r.p})
    res = pd.DataFrame(rows)
    (OUT / ds).mkdir(parents=True, exist_ok=True)
    res.to_csv(OUT / ds / "composition_null.tsv", sep="\t", index=False)
    if len(res):
        print(ds, f"composition test: {len(res)} null tests, share p < 0.05 = {(res.p < ALPHA).mean():.3f}")
    return 0


STATE_SIZES = (5, 10, 20, 40)


def run_state(ds, n_sets, seed, max_combos=4):
    """state_existence_test on random detected genes of one region x class: the test's p should be uniform (its null
    draws genes matched for expression and specificity) and the full rule (p < 0.05, tail ratio >= 2, donor spread
    >= 0.5) should pass rarely. Classes: those the project's state claims (contrast "within") name, most named first."""
    from oktopy import retest as rt
    rng = np.random.default_rng(seed)
    ctx = pd.read_csv(cfg.REPO / "claims" / "contexts.tsv", sep="\t").fillna("")
    named = ctx[ctx.contrast.eq("within")].cell_classes.str.split(";").explode().str.strip()
    present = set(rt.cell_obs(ds).cell_class.astype(str))
    classes = [c for c in named[named.ne("") & named.ne("ALL")].value_counts().index if c in present]
    regions = sorted(rt.cell_obs(ds).region.astype(str).unique())
    combos = [(r, c) for c in classes for r in regions][:max_combos]
    rows = []
    for reg, cls in combos:
        b = rt.class_cells(ds, (reg,), (cls,))
        if b is None:
            continue
        det = np.asarray((b.X > 0).mean(0)).ravel()
        pool = np.asarray(b.var_names)[det >= 0.01]
        for i in range(n_sets):
            size = int(STATE_SIZES[i % len(STATE_SIZES)])
            if len(pool) < size:
                continue
            r = rt.state_existence_test(ds, [reg], [cls], list(rng.choice(pool, size, replace=False)))
            if r is None:
                continue
            rows.append({"region": reg, "cell_class": cls, "size": size, "n_genes": r["n_genes"], "p": r["p"],
                         "z": r["z"], "tail_ratio": r["tail_ratio"], "donor_spread": r["donor_spread"],
                         "n_donors": r["n_case"], "passes": bool(r["p"] < ALPHA and r["tail_ratio"] >= 2
                                                                 and r["donor_spread"] >= 0.5)})
    res = pd.DataFrame(rows)
    (OUT / ds).mkdir(parents=True, exist_ok=True)
    res.to_csv(OUT / ds / "state_null.tsv", sep="\t", index=False)
    if len(res):
        print(ds, f"state existence on random genes: {len(res)} tests in {len(combos)} region x class, share p < 0.05 = "
              f"{(res.p < ALPHA).mean():.3f}, all rules passed = {res.passes.mean():.3f}")
    else:
        print(ds, "state existence: no region x class with enough unsorted nuclei (or no state claims): nothing tested")
    return 0


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / d
    return c - h, c + h


def status(share):
    return "calibrated" if BAND[0] <= share <= BAND[1] else ("liberal" if share > BAND[1] else "conservative")


SCORE_FILES = ["scores_null.tsv", "scores_small_null.tsv"]


def _read(d, names):
    parts = [pd.read_csv(d / n, sep="\t") for n in names if (d / n).is_file() and (d / n).stat().st_size > 5]
    return pd.concat(parts, ignore_index=True) if parts else None


def summary():
    rows = []
    for d in sorted(p for p in OUT.iterdir() if p.is_dir()) if OUT.is_dir() else []:
        f = d / "de_null.tsv"
        if f.is_file() and f.stat().st_size > 5:  # empty when the dataset has no DE fits
            t = pd.read_csv(f, sep="\t")
            if len(t):
                share = float(t.frac_p05.mean())
                rows.append({"dataset": d.name, "test": "edgeR genes", "n_null": len(t), "share_p05": share,
                             "ci_lo": np.nan, "ci_hi": np.nan, "median_lambda": float(t["lambda"].median()),
                             "fits_with_fdr05_hit": float((t.n_fdr05 > 0).mean()),
                             "median_n_units": float((t.n_case + t.n_ref).median()), "status": status(share)})
        # donor score rows are split by the genes actually tested (a 5-gene claim set can map to 4 genes)
        for name, files, sub in (
                ("competitive set, 2-4 genes (not primary)", ["sets_null.tsv"], lambda t: t[t.n_genes <= 4]),
                ("competitive set, 5+ genes (not primary)", ["sets_null.tsv"], lambda t: t[t.n_genes >= 5]),
                ("donor score, 2-4 genes", SCORE_FILES, lambda t: t[t.n_genes <= 4]),
                ("donor score, 5+ genes", SCORE_FILES, lambda t: t[t.n_genes >= 5]),
                ("composition CLR", ["composition_null.tsv"], lambda t: t),
                ("state existence, random genes", ["state_null.tsv"], lambda t: t),
                ("state existence, random genes, all rules", ["state_null.tsv"],
                 lambda t: t.assign(p=np.where(t.passes.astype(str).eq("True"), 0.0, 1.0)))):
            t = _read(d, files)
            if t is None:
                continue
            t = sub(t)
            if not len(t):
                continue
            k, n = int((t.p < ALPHA).sum()), len(t)
            lo, hi = wilson(k, n)
            units = (t.n_case + t.n_ref).median() if {"n_case", "n_ref"} <= set(t) else \
                t.n_donors.median() if "n_donors" in t else np.nan
            rows.append({"dataset": d.name, "test": name, "n_null": n, "share_p05": k / n, "ci_lo": lo, "ci_hi": hi,
                         "median_lambda": np.nan, "fits_with_fdr05_hit": np.nan, "median_n_units": units,
                         "status": status(k / n)})
    res = pd.DataFrame(rows)
    out = OUT.parent / "calibration_summary.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    res.round(4).to_csv(out, sep="\t", index=False)
    if len(res):
        print(res.round(3).to_string(index=False))
    print("wrote", out)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="oktopy calibrate", description=__doc__.split("\n\n")[0])
    ap.add_argument("what", choices=["de", "sets", "scores", "composition", "state", "summary"])
    ap.add_argument("datasets", nargs="*")
    ap.add_argument("--perms", type=int)
    ap.add_argument("--max-fits", type=int, default=40)
    ap.add_argument("--max-sets", type=int)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--small", action="store_true", help="scores: claims' gene sets of 2-4 genes instead of 5+")
    a = ap.parse_args(argv)
    if a.what == "summary":
        return summary()
    if not a.datasets:
        ap.error("name at least one dataset")
    codes = []
    for ds in a.datasets:
        (OUT / ds).mkdir(parents=True, exist_ok=True)
        if a.what == "de":
            codes.append(run_de(ds, a.perms or 5, a.max_fits, a.seed))
        elif a.what == "sets":
            codes.append(run_sets(ds, a.max_sets or 60, a.seed))
        elif a.what == "state":
            codes.append(run_state(ds, a.max_sets or 40, a.seed))
        elif a.what == "scores":
            codes.append(run_scores(ds, a.perms or 20, a.max_fits, a.max_sets or 30, a.seed, a.small))
        else:
            codes.append(run_composition(ds, a.perms or 100, a.seed))
    return max(codes)


if __name__ == "__main__":
    sys.exit(main())
