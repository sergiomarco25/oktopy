"""Power and informativeness, from the claim store (phase 7a; no verdict or store output changes).

  oktopy power [--db STORE] [--out DIR]

Writes to DIR (default results/robustness/):
  power_tests.tsv    one row per informative effect-family test (claim x dataset x variant): the power it had to
                     detect a true standardised effect of 0.2 (trivial), 0.5 (important) and the claim's pooled
                     estimate; `underpowered` = power at 0.5 below 0.5
  inconclusive.tsv   one row per Inconclusive claim: why (underpowered / heterogeneous / small_if_any /
                     pattern_not_detected), the pooled power at 0.5, and the donors per group a new two-group cohort
                     would need for 80% power at the pooled estimate
Effect sizes and standard errors follow claims/RATING.md section 2: for two groups SE = sqrt(1/n1 + 1/n2 +
d^2/(2(n1+n2))); for a correlation with n units, r = d/sqrt(d^2+4) and SE = sqrt(4/((n-1)(1-r^2))).
Power at alpha 0.05 two-sided, in the expected direction: Phi(|d|/SE - 1.96).
"""
import argparse
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from oktopy import config as cfg

Z = stats.norm.isf(0.025)
Z80 = stats.norm.ppf(0.8)
EFFECT_TEMPLATES = {"composition_clr", "composition_clr_sex", "composition_labels", "de_genes", "geneset_shift",
                    "state_genes_shift", "trajectory_genes", "communication_genes"}
TRIVIAL, IMPORTANT = 0.2, 0.5  # claims/RATING.md thresholds (equiv, important)
I2_SERIOUS = 0.5


def se_two_group(d, n1, n2):
    return np.sqrt(1 / n1 + 1 / n2 + d ** 2 / (2 * (n1 + n2)))


def se_correlation(d, n):
    r = d / np.sqrt(d ** 2 + 4)
    return np.sqrt(4 / ((n - 1) * (1 - r ** 2)))


def power(d, n1, n2=0):
    """Power to detect a true standardised effect d (two-sided alpha 0.05, expected direction)."""
    d = np.abs(np.asarray(d, float))
    n1, n2 = np.asarray(n1, float), np.asarray(n2, float)
    two = n2 > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        se = np.where(two, se_two_group(d, n1, np.where(two, n2, 1.0)), se_correlation(d, np.maximum(n1, 2.0)))
    return stats.norm.cdf(d / se - Z)


def donors_per_group(d, target=0.8):
    """Donors per group (two equal groups) for `target` power at effect d (normal approximation)."""
    d = np.abs(np.asarray(d, float))
    with np.errstate(divide="ignore"):
        return np.ceil(2 * (Z + stats.norm.ppf(target)) ** 2 / d ** 2)


def test_power(cd, claims):
    t = cd[cd.template_used.isin(EFFECT_TEMPLATES) & ~cd.outcome.eq("not_testable")].copy()
    t = t[t.n_case.notna()]
    t["n_ref"] = t.n_ref.fillna(0)
    t = t.merge(claims[["claim_id", "est", "verdict"]], on="claim_id", how="left")
    t["kind"] = np.where(t.n_ref > 0, "two_group", "correlation")
    t["power_02"] = power(TRIVIAL, t.n_case, t.n_ref)
    t["power_05"] = power(IMPORTANT, t.n_case, t.n_ref)
    t["power_est"] = np.where(t.est.notna(), power(t.est.fillna(0), t.n_case, t.n_ref), np.nan)
    t["underpowered"] = t.power_05 < 0.5
    cols = ["claim_id", "dataset", "role", "variant", "contrast", "template_used", "outcome", "n_case", "n_ref",
            "kind", "power_02", "power_05", "power_est", "underpowered"]
    return t[cols].round(4)


def inconclusive_reasons(claims):
    c = claims[claims.verdict.eq("Inconclusive")].copy()
    se = (c.hi - c.lo) / (2 * Z)
    c["pooled_power_05"] = stats.norm.cdf(IMPORTANT / se - Z)
    eff = c.family.eq("effect") & c.est.notna()
    het = eff & (c.k_cohorts.fillna(0) >= 2) & (c.i2.fillna(0) >= I2_SERIOUS)
    small = eff & ~het & (c.lo > -IMPORTANT) & (c.hi < IMPORTANT)
    c["reason"] = np.select([~eff, het, small], ["pattern_not_detected", "heterogeneous", "small_if_any"],
                            "underpowered")
    d = c.est.abs()
    c["donors_per_group_80"] = np.where(eff & (d >= 0.1), donors_per_group(d.where(d >= 0.1, 1)), np.nan)
    notes = {"heterogeneous": "datasets disagree (I2 >= 0.5): a single new cohort is unlikely to settle it",
             "small_if_any": "the interval excludes important effects in both directions: any effect is small",
             "underpowered": "the interval spans no effect and an important one: more donors would settle it",
             "pattern_not_detected": "pattern family: too few datasets detect it; power not modelled"}
    c["note"] = c.reason.map(notes)
    near_zero = eff & (d < 0.1)
    c.loc[near_zero, "note"] = c.loc[near_zero, "note"] + "; pooled estimate near zero"
    cols = ["claim_id", "origin", "family", "k_cohorts", "donors", "est", "lo", "hi", "i2", "reason",
            "pooled_power_05", "donors_per_group_80", "note"]
    return c[cols].round(4)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="oktopy power", description=__doc__.split("\n\n")[0])
    ap.add_argument("--db", default=str(cfg.STORE_DB))
    ap.add_argument("--out", default=str(cfg.RESULTS / "robustness"))
    a = ap.parse_args(argv)
    con = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
    claims = pd.read_sql("select claim_id, origin, verdict, family, k_cohorts, donors, est, lo, hi, i2 from claim", con)
    cd = pd.read_sql("select claim_id, dataset, role, variant, contrast, template_used, outcome, n_case, n_ref "
                     "from claim_dataset", con)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tp, inc = test_power(cd, claims), inconclusive_reasons(claims)
    tp.to_csv(out / "power_tests.tsv", sep="\t", index=False)
    inc.to_csv(out / "inconclusive.tsv", sep="\t", index=False)
    print(f"{len(tp)} effect-family tests: median power at 0.5 = {tp.power_05.median():.2f}; "
          f"{tp.underpowered.mean():.0%} underpowered (power at 0.5 < 0.5)")
    print(f"{len(inc)} Inconclusive claims by reason:")
    print(inc.groupby(["origin", "reason"]).size().to_string())
    print("wrote", out / "power_tests.tsv", "and", out / "inconclusive.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
