"""Claim re-tests switched on by curated tables: the share of author labels (claims/state_labels.tsv) and the
sex-stratified tests (claims/sex_claims.tsv, [tests] sex). Without the tables nothing changes."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ad = pytest.importorskip("anndata")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from oktopy import retest  # noqa: E402


def write_cells(tmp_path, monkeypatch, n_donors=12):
    """Donors d0..: the first half are cases; cases have 40% of their microglia in state S1, references 10%."""
    rows = []
    for i in range(n_donors):
        case = i < n_donors // 2
        n_s1 = 40 if case else 10
        rows += [{"donor": f"d{i}", "region": "R", "disease": "X" if case else "control", "sex": "F" if i % 2 else "M",
                  "cell_class": "Micro", "state": "S1" if j < n_s1 else "S0"} for j in range(100)]
        rows += [{"donor": f"d{i}", "region": "R", "disease": "X" if case else "control", "sex": "F" if i % 2 else "M",
                  "cell_class": "Astro", "state": "A0"} for _ in range(50)]
    obs = pd.DataFrame(rows)
    obs.index = obs.index.astype(str)
    ad.AnnData(np.zeros((len(obs), 1), np.float32), obs=obs).write_h5ad(tmp_path / "ds.h5ad")
    monkeypatch.setattr(retest, "H", tmp_path)
    monkeypatch.setattr(retest, "UNIT_COLUMNS", ["donor", "region", "disease", "sex", "cell_class"])
    monkeypatch.setattr(retest, "unit_groups", lambda meta: pd.DataFrame(
        {"c": np.where(meta.disease.astype(str).eq("X"), "case", "ref")}, index=meta.index))
    retest.label_obs.cache_clear()
    retest.cell_obs.cache_clear()


def test_no_tables_no_tests(tmp_path, monkeypatch):
    monkeypatch.setattr(retest, "STATE_LABELS_FILE", tmp_path / "none.tsv")
    monkeypatch.setattr(retest, "SEX_CLAIMS_FILE", tmp_path / "none.tsv")
    retest.state_labels.cache_clear()
    retest.sex_claims.cache_clear()
    assert retest.state_labels().empty and retest.sex_claims().empty


def test_label_share(tmp_path, monkeypatch):
    write_cells(tmp_path, monkeypatch)
    spec = pd.Series({"label_column": "state", "labels": "S1", "within": "Micro", "measure": ""})
    df = retest.label_composition_test("ds", ["R"], spec, "c")
    assert len(df) == 1 and df.z.iloc[0] > 0 and df.p.iloc[0] < 0.01
    assert df.prop_case.iloc[0] == pytest.approx(0.4) and df.prop_ref.iloc[0] == pytest.approx(0.1)
    res = retest.summarise(df, 1, "composition_labels")
    assert res["outcome"] == "reproduced"
    # the astrocytes are not in the denominator: with all nuclei the shares halve but the test still holds
    spec["within"] = ""
    df = retest.label_composition_test("ds", ["R"], spec, "c")
    assert df.prop_case.iloc[0] == pytest.approx(40 / 150)


def test_sex_of():
    assert list(retest.sex_of(["F", "male", "M", "Female", "unknown"])[:4]) == ["female", "male", "male", "female"]
    assert pd.isna(retest.sex_of(["unknown"])[0])


def test_pathology_shift(monkeypatch):
    rng = np.random.default_rng(1)
    units = [f"u{i}" for i in range(20)]
    stage = np.tile(np.arange(5), 4).astype(float)
    genes = pd.DataFrame(rng.normal(size=(200, 20)) + 0.8 * (stage - stage.mean()), columns=units)
    monkeypatch.setattr(retest, "unit_logcpm_z", lambda ds, reg, cls: genes)
    monkeypatch.setattr(retest, "unit_table", lambda ds: pd.DataFrame({"m": stage.astype(str)}, index=units))
    up = retest.pathology_shift_test("ds", ["R"], ["Oligo"], units, "m")
    assert up.z.iloc[0] > 3 and up.p.iloc[0] < 0.01 and up.effect.iloc[0] > 0
    flat = pd.DataFrame(rng.normal(size=(200, 20)), columns=units)
    monkeypatch.setattr(retest, "unit_logcpm_z", lambda ds, reg, cls: flat)
    assert retest.pathology_shift_test("ds", ["R"], ["Oligo"], units, "m").p.iloc[0] > 0.01
    assert retest.pathology_shift_test("ds", ["R"], ["Oligo"], units[:6], "m").empty  # too few units
    assert retest.pathology_shift_test("ds", ["R"], ["Oligo"], units, "absent").empty  # no such measure


def test_label_claim_elsewhere(tmp_path, monkeypatch):
    """A claim whose labels exist in one dataset: other datasets cannot test it when other_datasets is "none"."""
    f = tmp_path / "state_labels.tsv"
    f.write_text("claim_id\tdataset\tlabel_column\tlabels\twithin\tmeasure\tother_datasets\tnote\n"
                 "p-C1\tds\tstate\tS1\tMicro\t\tnone\tx\n")
    monkeypatch.setattr(retest, "STATE_LABELS_FILE", f)
    retest.state_labels.cache_clear()
    ctx = pd.Series({"claim_id": "p-C1", "regions": "R", "cell_classes": "Micro", "direction": "NA",
                     "contrast": "c", "test_template": "composition_donor", "cell_state": "S1"})
    res = retest._test_claim(ctx, "other", pd.DataFrame(columns=["gene", "kind"]))
    assert res[0]["outcome"] == "not_testable" and "only in ds" in res[0]["note"]
    retest.state_labels.cache_clear()



def test_both_sexes_one_row():
    """A claim about both sexes gives one row per dataset, whatever the row order: the sex that tests it hardest."""
    f = {"outcome": "reproduced", "z": -1.0, "p": 0.32, "note": None}
    m = {"outcome": "contradicted", "z": -2.5, "p": 0.012}
    na = {"outcome": "not_testable", "note": "too few units"}
    for rows in ({"female": f, "male": m}, {"male": m, "female": f}):
        r = retest.both_sexes(rows, null=True)  # no-difference claim: the larger change, p doubled
        assert r["variant"] == "both" and r["other_sex_z"] == -1.0 and r["p"] == pytest.approx(0.024)
        assert r["z"] == pytest.approx(-2.257, abs=1e-3) and r["note"].startswith("both sexes: male")
    assert retest.both_sexes({"female": f, "male": m}, null=False)["other_sex_p"] == 0.32  # weaker support: male
    one = retest.both_sexes({"female": na, "male": f}, null=True)  # one sex testable: its own p
    assert one["p"] == 0.32 and one["variant"] == "both" and one["note"] == "both sexes: male shown"
    assert retest.both_sexes({"female": na, "male": na}, null=True)["outcome"] == "not_testable"
