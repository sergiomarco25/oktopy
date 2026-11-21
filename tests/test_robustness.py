"""Phase 7a: method registry, shuffled-label calibration script, power.

  python3 -I tests/test_robustness.py   or   pytest tests/
The registry and R checks need only the standard library; the power checks need numpy, pandas and scipy.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from oktopy import registry  # noqa: E402

ENGINE = ROOT / "src" / "oktopy"


def test_registry_is_valid_for_this_project():
    reg = registry.load()
    assert registry.validate(reg, registry.templates()) == []


def test_registry_rules_are_enforced():
    reg = registry.load()
    task = reg["task"]["de.case_control"]
    task["method"].append(dict(task["method"][0], id="second_primary"))
    errs = registry.validate(reg, registry.templates())
    assert any("primaries" in e for e in errs)
    reg = registry.load()
    reg["task"]["de.case_control"]["method"][0]["disabled"] = True
    assert any("cannot switch off the primary" in e for e in registry.validate(reg))
    reg = registry.load()
    reg["task"]["de.case_control"]["method"][1]["tier"] = "favourite"
    assert any("tier must be" in e for e in registry.validate(reg))


def test_every_primary_cites_the_book_or_a_benchmark_and_says_why():
    for r in registry.rows(registry.load()):
        assert r["evidence"] and r["why"], r["id"]


def test_markdown_lists_every_task():
    reg = registry.load()
    md = registry.markdown(reg)
    for name in reg["task"]:
        assert f"(`{name}`)" in md


def _fit_one(path):
    s = path.read_text()
    a = s.index("fit_one <- function(")
    b = s.index("\n}\n", a)
    return s[a:b + 2]


def test_calibration_refits_with_the_de_cache_model():
    assert _fit_one(ENGINE / "de_calibration.R") == _fit_one(ENGINE / "de_cache.R")


def _covariate_rule(path):
    s = path.read_text()
    a = s.index("  covs <- c()\n")
    b = s.index("    covs <- head(covs, -1)\n  }\n", a)
    return s[a:b]


def test_sensitivity_methods_use_the_de_cache_design():
    """Phase 7b: the sensitivity DE methods build the DE cache's design (same covariates, same gene filter)."""
    assert _covariate_rule(ENGINE / "de_sensitivity.R") == _covariate_rule(ENGINE / "de_cache.R")
    s = (ENGINE / "de_sensitivity.R").read_text()
    assert "keep <- filterByExpr(y, design)" in s and "calcNormFactors(y[keep, , keep.lib.sizes = FALSE])" in s


def test_multiverse_set_keys_and_rows():
    import pandas as pd
    from oktopy import multiverse as mv
    assert mv.set_key({"A": 1, "B": -1}) == mv.set_key({"B": -1, "A": 1})
    assert mv.set_key({"A": 0, "B": -1}) == mv.set_key({"A": 1, "B": -1})      # no direction is tested as up
    assert mv.set_key({"A": 1, "B": 1}) != mv.set_key({"A": 1, "B": -1})
    df = pd.DataFrame({"set_test": ["single gene (FDR)", "donor signature score (Mann-Whitney)", "", ""],
                       "template_used": ["de_genes", "state_genes_shift", "composition_clr", "state_existence"]})
    assert mv.applies(df, "deseq2").tolist() == [True, False, False, False]
    assert mv.applies(df, "camera").tolist() == [False, True, False, False]
    assert mv.applies(df, "propeller").tolist() == [False, False, True, False]
    assert not mv.method_ran("competitive set test (random sets)", "camera")      # a fallback is not the method
    assert mv.method_ran(mv.LABELS["camera"], "camera")


def test_adjusted_score_recovers_a_shift_with_covariates():
    import numpy as np
    import pandas as pd
    from oktopy import multiverse as mv
    rng = np.random.default_rng(0)
    n = 40
    u = pd.DataFrame({"age": rng.normal(80, 8, n).round(), "sex": rng.choice(["M", "F"], n), "pmi": [""] * n})
    x = np.r_[np.ones(20), np.zeros(20)]
    y = 0.8 * x + 0.05 * (u.age.values - 80) + rng.normal(0, 0.5, n)
    beta, p, covs = mv.ols_term(y, x, u)
    assert covs == ["age_s", "sex"] and 0.5 < beta < 1.1 and p < 0.001   # PMI unknown for all: not a covariate
    beta0, p0, _ = mv.ols_term(rng.normal(0, 1, n), x, u)
    assert p0 > 0.001


def test_adjusted_score_native_effect_matches_the_rating_scale():
    """Without covariates the adjusted test's own effect size equals what rating.to_smd makes of an unadjusted test;
    with a covariate collinear with the group, the standard error grows and the effect is not shrunk. A test that
    reports no effect size of its own leaves the re-test row as it was."""
    import numpy as np
    import pandas as pd
    from oktopy import multiverse as mv
    from oktopy import rating
    from oktopy import retest as rt
    rng = np.random.default_rng(2)
    n = 60
    none = pd.DataFrame({"age": [""] * n, "sex": [""] * n, "pmi": [""] * n})  # no covariate known for all units
    x = np.r_[np.ones(30), np.zeros(30)]
    y = 0.6 * x + rng.normal(0, 1, n)
    f = mv.ols_fit(y, x, none)
    d, se = mv.native_smd(f, x, stage=False)
    pooled = np.sqrt(((y[x == 1].var(ddof=1) + y[x == 0].var(ddof=1)) / 2))
    assert abs(d - (y[x == 1].mean() - y[x == 0].mean()) / pooled) < 1e-9         # Cohen's d
    d_z, se_z, _ = rating.to_smd(f["beta"] / f["se"], 30, 30, 1)                  # t as z: close, not equal
    assert abs(d - d_z) < 0.02 and abs(se - se_z) < 0.01
    s = np.repeat(np.arange(6.0), 10)
    ys = 0.3 * s + rng.normal(0, 1, n)
    d, se = mv.native_smd(mv.ols_fit(ys, s, none), s, stage=True)
    r = np.corrcoef(ys, s)[0, 1]
    d_r, se_r, _ = rating.to_smd(r * np.sqrt(n), n, 0, 1)
    assert abs(d - d_r) < 0.05 * abs(d_r) and abs(se - se_r) < 0.05 * se_r
    age = 75 + 8 * (x - 0.5) * 2 + rng.normal(0, 3, n)                             # age tracks the group (r ~ 0.9)
    u = pd.DataFrame({"age": age.round(), "sex": [""] * n, "pmi": [""] * n})
    ya = 0.6 * x + rng.normal(0, 1, n)
    fa, f0 = mv.ols_fit(ya, x, u), mv.ols_fit(ya, x, none)
    da, sea = mv.native_smd(fa, x, stage=False)
    d0, se0 = mv.native_smd(f0, x, stage=False)
    assert fa["covs"] == ["age_s"] and sea > 1.5 * se0
    d_conv = rating.to_smd(fa["beta"] / fa["se"], 30, 30, 1)[0]                   # z converted: shrunk by ~sqrt(VIF)
    assert abs(da) > 1.5 * abs(d_conv)
    df = pd.DataFrame({"z": [2.0, 1.0], "p": [0.05, 0.3], "effect": [0.5, 0.2], "n_case": [10, 10],
                       "n_ref": [10, 10], "region": ["A", "B"], "cell_class": ["C", "C"], "test": ["t", "t"]})
    assert "d_native" not in rt.summarise(df, 1, "x")
    out = rt.summarise(df.assign(d_native=[0.8, 0.4], se_native=[0.5, 0.3]), 1, "x")
    assert abs(out["d_native"] - 0.6) < 1e-9 and abs(out["se_native"] - 0.4) < 1e-9


def _rscript_with(*pkgs):
    """The engine's Rscript if it loads every package, else None (the R checks are then skipped)."""
    import subprocess
    from oktopy import de_cache
    r = de_cache.rscript()
    code = ";".join(f"suppressPackageStartupMessages(library({p}))" for p in pkgs)
    try:
        return r if subprocess.run([r, "--vanilla", "-e", code], capture_output=True).returncode == 0 else None
    except OSError:
        return None


def test_sensitivity_de_methods_fit_group_and_stage():
    """Every DE sensitivity method writes a group and a stage table on a small synthetic DE-cache input; the Wilcoxon
    method's stage test is Spearman's rho (a transposed rank matrix once made every stage fit fail)."""
    import gzip
    import subprocess
    import tempfile
    import numpy as np
    import pandas as pd
    from scipy import stats
    r = _rscript_with("data.table", "edgeR", "limma")
    if r is None:
        print("skipped: no Rscript with data.table, edgeR and limma")
        return
    methods = ["edger_ql", "limma_voom", "wilcoxon"] + (["deseq2"] if _rscript_with("DESeq2") else [])
    rng = np.random.default_rng(1)
    n, genes = 12, [f"G{i}" for i in range(300)]
    stage = np.array([0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 6, 6])
    units = pd.DataFrame({"unit": [f"d{i}|R" for i in range(n)], "age": rng.integers(60, 95, n),
                          "sex": ["M", "F"] * 6, "pmi": rng.uniform(2, 30, n).round(1), "stage": stage,
                          "AD_vs_control": np.where(stage >= 3, "case", "ref")})
    lam = rng.uniform(20, 200, len(genes))[:, None] * np.ones(n)
    lam[0] *= np.exp(0.4 * stage)                                         # G0 rises with stage
    counts = pd.DataFrame(rng.poisson(lam), index=genes, columns=units.unit)
    man = pd.DataFrame([{"dataset": "syn", "region": "R", "cell_class": "C", "file": "R__C.csv.gz",
                         "contrast": "AD_vs_control", "n_case": 6, "n_ref": 6, "median_cells": 100, "run": True}])
    with tempfile.TemporaryDirectory() as tmp:
        inp = Path(tmp) / "syn" / "input"
        inp.mkdir(parents=True)
        units.to_csv(inp / "units.tsv", sep="\t", index=False)
        man.to_csv(inp / "manifest.tsv", sep="\t", index=False)
        with gzip.open(inp / "R__C.csv.gz", "wt") as fh:
            counts.to_csv(fh)
        for m in methods:
            out = Path(tmp) / "out" / m
            out.mkdir(parents=True)
            subprocess.run([r, "--vanilla", str(ENGINE / "de_sensitivity.R"), m, str(inp.parent), str(out), "-",
                            "fit", "0", "0", "0", "ad_pathology", "stage", "", ""], check=True, capture_output=True)
            for con in ("AD_vs_control", "ad_pathology"):
                t = pd.read_csv(out / f"de_{m}" / f"R__C__{con}.tsv.gz", sep="\t").set_index("gene")
                assert t.p.between(0, 1).all() and len(t) >= 50, (m, con)
                assert t.loc["G0", "logFC"] > 0 and t.loc["G0", "p"] < (0.05 if con == "AD_vs_control" else 1e-4), \
                    (m, con)                       # the group split leaves the stage gradient inside each group
            if m == "wilcoxon":
                cpm = counts / counts.sum(0) * 1e6
                rho = stats.spearmanr(cpm.loc["G1"], stage)[0]
                t = pd.read_csv(out / "de_wilcoxon" / "R__C__ad_pathology.tsv.gz", sep="\t").set_index("gene")
                assert abs(t.loc["G1", "logFC"] - rho) < 0.15                  # TMM vs plain CPM: close, not equal


def test_calibration_shuffles_labels_not_covariates():
    s = (ENGINE / "de_calibration.R").read_text()
    assert re.search(r"uu\$group <- factor\(sample\(g\[sel\]\)", s)
    assert "set.seed(seed)" in s


def test_power_matches_the_rating_effect_model():
    import numpy as np
    from oktopy import power
    # 10 vs 10 donors: SE = sqrt(0.2 + 0.25/40); power = Phi(0.5/SE - 1.96)
    assert abs(float(power.power(0.5, 10, 10)) - 0.195) < 0.005
    assert abs(float(power.power(0.5, 50, 50)) - 0.69) < 0.01
    assert abs(float(power.power(0.5, 100, 100)) - 0.936) < 0.005
    assert float(power.power(0.5, 100, 100)) > float(power.power(0.2, 100, 100))
    assert 0 < float(power.power(0.5, 30, 0)) < 1                      # correlation with 30 units
    assert int(power.donors_per_group(0.5)) == 63
    assert np.isnan(power.donors_per_group(np.nan))


def test_inconclusive_reasons():
    import pandas as pd
    from oktopy import power
    c = pd.DataFrame([
        {"claim_id": "a", "origin": "published", "verdict": "Inconclusive", "family": "effect", "k_cohorts": 1,
         "donors": 20, "est": 0.3, "lo": -0.6, "hi": 1.2, "i2": 0.0},
        {"claim_id": "b", "origin": "published", "verdict": "Inconclusive", "family": "effect", "k_cohorts": 4,
         "donors": 300, "est": 0.3, "lo": -0.4, "hi": 1.0, "i2": 0.7},
        {"claim_id": "c", "origin": "published", "verdict": "Inconclusive", "family": "effect", "k_cohorts": 5,
         "donors": 600, "est": 0.15, "lo": -0.05, "hi": 0.35, "i2": 0.1},
        {"claim_id": "d", "origin": "published", "verdict": "Inconclusive", "family": "pattern", "k_cohorts": 3,
         "donors": 90, "est": None, "lo": None, "hi": None, "i2": None},
        {"claim_id": "e", "origin": "published", "verdict": "Supported", "family": "effect", "k_cohorts": 3,
         "donors": 90, "est": 0.6, "lo": 0.3, "hi": 0.9, "i2": 0.0},
    ])
    r = power.inconclusive_reasons(c).set_index("claim_id").reason.to_dict()
    assert r == {"a": "underpowered", "b": "heterogeneous", "c": "small_if_any", "d": "pattern_not_detected"}


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
