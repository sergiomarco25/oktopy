"""Harmonisation stops when nuclei lack a donor or region: pooled, they would form one false donor or region."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ad = pytest.importorskip("anndata")
pytest.importorskip("scanpy")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from oktopy import harmonise  # noqa: E402


@pytest.mark.parametrize("key", ["donor", "region"])
def test_missing_unit_key_stops(tmp_path, monkeypatch, key):
    monkeypatch.setattr(harmonise, "OUT", tmp_path)
    obs = pd.DataFrame({"donor": ["d1", "d2", "d3"], "region": ["R"] * 3}, index=["a", "b", "c"])
    obs.loc["c", key] = np.nan
    a = ad.AnnData(np.ones((3, 2), np.float32), obs=obs)
    with pytest.raises(AssertionError, match=f"1 nuclei without a {key}"):
        harmonise.finalise(a, "ds", {})
    assert not (tmp_path / "ds.h5ad").exists()
