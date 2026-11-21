"""Linked-data exports of a release (oktopy.linked_data), on the example store.

Standard library only:  python3 -I tests/test_linked_data.py   or   pytest tests/
"""
import json
import shutil
import sqlite3
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCHEMA_ORG = {"about", "alternateName", "appearance", "bestRating", "claimReviewed", "datePublished", "identifier",
              "inDefinedTermSet", "isBasedOn", "isPartOf", "itemReviewed", "keywords", "name", "ratingExplanation",
              "ratingValue", "reviewAspect", "reviewBody", "reviewRating", "termCode", "text", "worstRating"}


def _module():
    from oktopy import linked_data
    return linked_data


def _terms(x, out):
    if isinstance(x, dict):
        for k, v in x.items():
            if not k.startswith("@"):
                out.add(k)
            _terms(v, out)
    elif isinstance(x, list):
        for v in x:
            _terms(v, out)
    return out


def test_claims_jsonld_covers_every_claim(store_db):
    doc = _module().claims_jsonld(store_db)
    n = sqlite3.connect(f"file:{store_db}?mode=ro", uri=True).execute("SELECT count(*) FROM claim").fetchone()[0]
    assert doc["@context"] == {"@vocab": "https://schema.org/"}
    assert len(doc["@graph"]) == n
    assert _terms(doc, set()) <= SCHEMA_ORG  # only the schema.org properties we checked by hand
    r = doc["@graph"][0]
    assert r["@type"] == "ClaimReview" and r["itemReviewed"]["@type"] == "Claim"
    # the example's fictitious papers have no DOI and its cell types no ontology ids: no links are made up
    assert all(a.get("@id", "").startswith("https://doi.org/") for a in [r["itemReviewed"].get("appearance") or {}] if a)


def test_ro_crate_lists_the_release(store_db, tmp_path=None):
    ld = _module()
    d = Path(tmp_path or tempfile.mkdtemp())
    shutil.copy(store_db, d / "claim_store.sqlite")
    (d / "tables").mkdir()
    (d / "tables" / "claim.tsv").write_text("claim_id\n")
    ld.write(d, {"name": "test", "version": "0.0.1", "date": "2026-01-01"})
    crate = json.loads((d / "ro-crate-metadata.json").read_text())
    ids = {n["@id"]: n for n in crate["@graph"]}
    assert ids["ro-crate-metadata.json"]["about"] == {"@id": "./"}
    root = ids["./"]
    assert root["license"] == {"@id": ld.CC_BY} and root["author"]
    assert {"claim_store.sqlite", "claims.jsonld", "tables/", "release.json", "SHA256SUMS"} <= {p["@id"] for p in root["hasPart"]}
    assert ids["tables/"]["hasPart"] == [{"@id": "tables/claim.tsv"}]
    assert ids["claim_store.sqlite"]["encodingFormat"] == "application/vnd.sqlite3"
