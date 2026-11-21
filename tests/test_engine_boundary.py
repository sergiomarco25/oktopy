"""The engine (src/oktopy) is field-independent: it names no disease, tissue, dataset or paper and imports no
project package. Facts about a field come from its project's config.

Standard library only:  pytest tests/test_engine_boundary.py
"""
import re
from pathlib import Path

from conftest import EXAMPLE_PAPERS

ROOT = Path(__file__).resolve().parents[1]
ENGINE_DIR = ROOT / "src" / "oktopy"
# the engine's modules; example.py is the worked example itself (it names its own papers)
ENGINE = sorted(p for p in ENGINE_DIR.glob("*.py") if p.name != "example.py")
ENGINE_FILES = ENGINE + sorted(ENGINE_DIR.glob("*.R")) + sorted(ENGINE_DIR.glob("*.html")) + sorted(
    ENGINE_DIR.glob("*.toml"))
# words that belong to one field, never to the engine (the first two fields: Parkinson's and Alzheimer's disease);
# case-sensitive, so the pandas alias pd does not count
FIELD = (r"Parkinson|\bPD\b|\bPDD\b|\bDLB\b|Lewy|\blewy|nigra|dopamin|Braak|\bbraak|striat|putamen|caudate|"
         r"substantia|synuclein|\bSNCA\b|LRRK2|\bGBA\b|Alzheimer|\bAPOE\b|amyloid|\bMSNs?\b")


def test_engine_names_nothing_of_a_field():
    pat = re.compile(FIELD + "|" + "|".join(rf"\b{re.escape(i)}\b" for i in sorted(EXAMPLE_PAPERS)))
    hits = [f"{p.name}:{i}: {m.group(0)}" for p in ENGINE_FILES
            for i, line in enumerate(p.read_text().splitlines(), 1) for m in [pat.search(line)] if m]
    assert not hits, f"field-specific names in the engine (move them to the project config): {hits}"


def test_engine_imports_no_project_package():
    for p in ENGINE:
        src = p.read_text()
        bad = re.findall(r"^\s*(?:from|import)\s+(\w+)", src, re.M)
        assert not {"claimstore", "example_claims"} & set(bad), p.name
    assert len(ENGINE) >= 30, [p.name for p in ENGINE]
