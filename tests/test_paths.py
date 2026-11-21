"""Static checks that keep the code portable. Standard library only.

- no cluster paths or personal e-mail addresses in any file of the repository (paths belong in a project's
  config/project.toml or its environment variables);
- every engine file that uses `cfg.` imports the config before its first use.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SKIP = {".git", "__pycache__", ".pytest_cache", "build", "dist"}
FILES = [p for p in ROOT.rglob("*") if p.is_file() and not SKIP & set(p.relative_to(ROOT).parts)
         and not p.name.endswith((".egg-info", ".pyc")) and p.name != "test_paths.py"]
CLUSTER = re.compile(r"/(lustre|ictstr01|gpfs|scratch/[a-z]+\.[a-z]+)/")
EMAIL = re.compile(r"[\w.+-]+@[A-Za-z][\w-]*\.[A-Za-z][\w.-]*")  # not package@1.2.3 pins
EMAIL_OK = {"noreply@anthropic.com"}


def _text(p):
    try:
        return p.read_text()
    except (UnicodeDecodeError, OSError):
        return ""


def test_no_cluster_paths():
    hits = [f"{p.relative_to(ROOT)}:{i}" for p in FILES for i, line in enumerate(_text(p).splitlines(), 1)
            if CLUSTER.search(line)]
    assert not hits, f"cluster paths: {hits}"


def test_no_personal_email():
    hits = [f"{p.relative_to(ROOT)}:{i}: {m}" for p in FILES for i, line in enumerate(_text(p).splitlines(), 1)
            for m in EMAIL.findall(line) if m not in EMAIL_OK and not m.endswith(("example.org", "example.com"))
            and "@v" not in m and "git@github.com" not in m]
    assert not hits, f"e-mail addresses: {hits}"


def test_config_imported_before_use():
    bad = []
    for p in sorted((ROOT / "src" / "oktopy").glob("*.py")):
        lines = _text(p).splitlines()
        first_use = next((i for i, ln in enumerate(lines) if re.search(r"\bcfg\.", ln)), None)
        if first_use is None:
            continue
        imp = next((i for i, ln in enumerate(lines) if "import config as cfg" in ln), None)
        if imp is None or imp > first_use:
            bad.append(p.name)
    assert not bad, f"cfg used before `from oktopy import config as cfg`: {bad}"
