"""`oktopy new` makes a working project for another field in one command, and `oktopy example` the worked example.

Generates projects into temporary folders and checks that they hold nothing of another field (no disease terms,
cluster paths or the example's ids in a new project), that every placeholder is filled, that their own checks pass,
and that a new project with claims alone builds, rates, serves and shows a store. Needs pytest, pyyaml, jsonschema,
pandas and scipy.
"""
import json
import os
import re
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest

from oktopy import __version__, example, newproject

pytest.importorskip("jsonschema")
pytest.importorskip("yaml")
pytest.importorskip("pandas")

FOREIGN = re.compile(r"Parkinson|\bPD\b|Alzheimer|\bAD\b|Lewy|nigra|Braak|striat|synuclein|/" + "lustre/|/" + "ictstr01/|"
                     r"claimstore|smith2024|jones2025|\{\{\w+\}\}")


def _new(base, name="new-field", **kw):
    values = newproject.values_for(name, "New field claim store", "test project", **kw)
    return newproject.generate(base / name, values)


def _files(out):
    return [p for p in out.rglob("*") if p.is_file() and "__pycache__" not in p.parts]


def test_new_project_is_clean_and_passes_its_checks(tmp_path):
    out = _new(tmp_path)
    for rel in ("config/project.toml", "config/schema.toml", "config/contrasts.toml", "config/harmonise.toml",
                "config/contexts.toml", "claims/example2024.yaml", "claims/schema/claim.schema.json",
                "claims/RATING.md", "src/new_field/build/example2024.py", "tests/test_project.py",
                "workflow/Snakefile", "workflow/profiles/default/config.yaml", "workflow/profiles/slurm/config.yaml",
                ".github/workflows/ci.yml", ".gitignore", "requirements.txt", "LICENSE", "LICENSE-DATA.md",
                "CITATION.cff", "README.md", ".claude/skills/claim-replication/SKILL.md",
                ".claude/skills/claim-replication/references/extraction.md"):
        assert (out / rel).is_file(), rel
    assert not (out / "src" / "oktopy").exists() and not (out / "src" / "claimstore").exists(), "engine copied"
    hits = [f"{p.relative_to(out)}:{i}: {m.group(0)}" for p in _files(out) if p.suffix not in (".json",)
            for i, line in enumerate(p.read_text().splitlines(), 1) for m in [FOREIGN.search(line)] if m
            and not (p.name == "Snakefile" and m.group(0).startswith("{{"))]
    assert not hits, hits[:20]
    conf = tomllib.loads((out / "config" / "project.toml").read_text())
    assert conf["project"]["oktopy"] == __version__ and conf["project"]["env_prefix"] == "NEW_FIELD"
    assert newproject.REQUIREMENT in (out / "requirements.txt").read_text()
    assert newproject.run_checks(out) is True


def test_bad_names_are_refused(tmp_path):
    for kw in ({"name": "Bad Name"}, {"name": "ok-name", "package": "oktopy"}):
        with pytest.raises(SystemExit):
            newproject.values_for(kw.pop("name"), "t", **kw)


def test_existing_folder_is_not_replaced_unless_it_is_a_project(tmp_path):
    (tmp_path / "other").mkdir()
    (tmp_path / "other" / "notes.txt").write_text("keep me")
    values = newproject.values_for("x-field", "X")
    with pytest.raises(SystemExit):
        newproject.generate(tmp_path / "other", values)
    with pytest.raises(SystemExit):
        newproject.generate(tmp_path / "other", values, force=True)
    assert (tmp_path / "other" / "notes.txt").read_text() == "keep me"


def test_new_project_builds_rates_and_serves_a_store(tmp_path):
    """With the engine alone, a new project merges its claims, maps contexts, builds and rates a store, serves it
    and makes its atlas page."""
    pytest.importorskip("scipy")
    out = _new(tmp_path)
    env = {**os.environ, "OKTOPY_PROJECT": str(out), "NEW_FIELD_WORK": str(tmp_path / "work")}
    cli = [sys.executable, "-m", "oktopy"]
    for step in ("merge", "contexts", "store", "rating", "store", "atlas"):
        r = subprocess.run(cli + [step], cwd=out, env=env, capture_output=True, text=True)
        assert r.returncode == 0, (step, r.stdout[-800:], r.stderr[-800:])
    db = tmp_path / "work" / "results" / "store" / "claim_store.sqlite"
    msgs = [{"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {"protocolVersion": "2025-06-18"}},
            {"jsonrpc": "2.0", "id": 2, "method": "tools/call", "params": {"name": "overview", "arguments": {}}},
            {"jsonrpc": "2.0", "id": 3, "method": "tools/call",
             "params": {"name": "get_claim", "arguments": {"claim_id": "example2024-C01"}}}]
    r = subprocess.run(cli + ["serve"], cwd=out, env={**env, "NEW_FIELD_STORE_DB": str(db)},
                       input="\n".join(json.dumps(m) for m in msgs) + "\n", capture_output=True, text=True, timeout=120)
    replies = {x["id"]: x for x in map(json.loads, r.stdout.splitlines())}
    assert replies[1]["result"]["serverInfo"]["name"] == "new-field"
    for i in (2, 3):
        assert not replies[i]["result"]["isError"], replies[i]["result"]["content"][0]["text"][:300]
    claim = json.loads(replies[3]["result"]["content"][0]["text"])["claim"]
    assert claim["verdict"] == "Untested"
    page = (out / "reports" / "atlas.html").read_text()
    assert "<title>New field claim store</title>" in page and "__TITLE__" not in page and "/*__DATA__*/null" not in page
    data = json.loads(page.split("const D = ", 1)[1].split(";\n", 1)[0])
    assert [c["id"] for c in data["claims"]] == ["example2024-C01"] and data["field"]["title"] == "New field claim store"


def test_project_finds_its_folder_from_a_subfolder(tmp_path):
    out = _new(tmp_path)
    sub = out / "claims"
    env = {k: v for k, v in os.environ.items() if k not in ("OKTOPY_PROJECT", "OKTOPY_CONFIG", "CLAIMSTORE_CONFIG")}
    r = subprocess.run([sys.executable, "-m", "oktopy", "config", "--show"], cwd=sub, env=env, capture_output=True,
                       text=True)
    assert r.returncode == 0 and f"project     {out}" in r.stdout, r.stdout + r.stderr


def test_config_json_for_workflow_tools(tmp_path):
    """The Snakefile reads paths and datasets through `oktopy config --json`, so Snakemake needs no oktopy import."""
    out = example.make(tmp_path / "ex")
    r = subprocess.run([sys.executable, "-m", "oktopy", "--project", str(out), "config", "--json"],
                       capture_output=True, text=True)
    c = json.loads(r.stdout)
    assert c["REPO"] == str(out) and c["DATASETS"] == ["smith2024", "jones2025"]
    assert c["HARMONISED"] == str(out / "data" / "harmonised") and c["MEM_GB"] == {"smith2024": 4, "jones2025": 4}
    assert "config --json" in (out / "workflow" / "Snakefile").read_text()


def test_example_project_passes_its_checks(tmp_path):
    out = example.make(tmp_path / "example")
    assert not (out / "claims" / "example2024.yaml").exists()
    assert {p.stem for p in (out / "claims").glob("*20*.yaml")} == {"smith2024", "jones2025"}
    conf = tomllib.loads((out / "config" / "project.toml").read_text())
    assert conf["pipeline"]["datasets"] == ["smith2024", "jones2025"]
    assert newproject.run_checks(out) is True


def test_example_shares_the_template_claim_vocabularies(tmp_path):
    """Only the field's regions and cell types differ, so claims stay comparable across fields."""
    tpl = tomllib.loads((newproject.TEMPLATE / "config" / "schema.toml").read_text())
    ex = tomllib.loads((example.make(tmp_path / "ex") / "config" / "schema.toml").read_text())
    for k in tpl["vocab"]:
        if k not in ("region", "cell_type"):
            assert ex["vocab"][k] == tpl["vocab"][k], k
    assert ex["required"] == tpl["required"] and ex["text"] == tpl["text"]
