"""The worked example end to end: `oktopy example` then `oktopy run` (harmonise, pseudobulk, edgeR DE, re-tests,
meta-analysis, store, rating, atlas), and the verdicts of the claims whose answer the synthetic data plant.

Needs oktopy[analysis] and R with edgeR and data.table; skipped without them unless OKTOPY_E2E=1 (CI sets it, so a
missing dependency fails there instead of skipping). Takes a few minutes.
"""
import os
import shutil
import sqlite3
import subprocess
import sys

import pytest

from oktopy import example

REQUIRED = os.environ.get("OKTOPY_E2E") == "1"

# claim -> verdicts the planted effects allow (synthetic data: _synthetic.py in the example's loaders)
EXPECTED = {
    "smith2024-C01": {"Supported"},                     # immune share x2 in cases
    "jones2025-C01": {"Supported"},
    "smith2024-C02": {"Supported"},                     # SPP1 x3 in immune cells
    "smith2024-C03": {"Supported"},                     # interferon genes x2 / x1.6 in epithelial cells
    "jones2025-C03": {"Supported"},
    "smith2024-C04": {"Contradicted"},                  # HIF1A claimed up, planted down
    "smith2024-C06": {"Supported"},                     # MKI67 rises with stage
    "jones2025-C02": {"Inconclusive", "No effect"},     # COL3A1 claimed down, unchanged
    "smith2024-C05": {"Inconclusive", "Supported (no effect)"},  # endothelial share unchanged (needs power)
    "jones2025-C04": {"Untested"},                      # spatial claim, no spatial data
}


def _have_r():
    r = shutil.which("Rscript")
    if r is None:
        return False
    code = "suppressPackageStartupMessages({library(edgeR); library(data.table)})"
    return subprocess.run([r, "--vanilla", "-e", code], capture_output=True).returncode == 0


def _skip_or_fail(reason):
    if REQUIRED:
        pytest.fail(reason)
    pytest.skip(reason)


def test_worked_example_runs_and_rates_the_planted_claims(tmp_path):
    for mod in ("anndata", "scanpy"):
        try:
            __import__(mod)
        except ImportError:
            _skip_or_fail(f"{mod} missing (pip install 'oktopy[analysis]')")
    if not _have_r():
        _skip_or_fail("no Rscript with edgeR and data.table")
    out = example.make(tmp_path / "example")
    env = {**os.environ, "OKTOPY_PROJECT": str(out)}
    r = subprocess.run([sys.executable, "-m", "oktopy", "run"], cwd=out, env=env, capture_output=True, text=True,
                       timeout=3600)
    assert r.returncode == 0, r.stdout[-4000:] + r.stderr[-4000:]
    db = out / "data" / "results" / "store" / "claim_store.sqlite"
    con = sqlite3.connect(db)
    got = dict(con.execute("SELECT claim_id, verdict FROM claim WHERE origin = 'published'"))
    assert set(got) == set(EXPECTED)
    wrong = {c: (got[c], sorted(ok)) for c, ok in EXPECTED.items() if got[c] not in ok}
    assert not wrong, wrong
    n_tests = con.execute("SELECT count(DISTINCT dataset) FROM claim_dataset WHERE outcome != 'not_testable'").fetchone()[0]
    assert n_tests == 2
    assert (out / "reports" / "atlas.html").stat().st_size > 50_000
    # certainty: two cohorts start at Moderate; never above High, the untested claim Very low
    cert = dict(con.execute("SELECT claim_id, certainty FROM claim WHERE origin = 'published'"))
    assert cert["jones2025-C04"] == "Very low" and {cert[c] for c in got} <= {"High", "Moderate", "Low", "Very low"}
