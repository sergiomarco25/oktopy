"""What a fresh clone and an installed package hold: every file of the project template, the worked example and the
skill is tracked by git (an ignore rule once left the example's loaders out of the repository), and the package-data
globs of pyproject.toml cover them. Standard library only; skipped outside a git checkout."""
import fnmatch
import shutil
import subprocess
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "src" / "oktopy"
DATA_DIRS = ("template", "example", "skill")


def _data_files():
    return sorted(p.relative_to(ROOT).as_posix() for d in DATA_DIRS for p in (PKG / d).rglob("*")
                  if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc")


def test_package_data_is_tracked_by_git():
    if shutil.which("git") is None or not (ROOT / ".git").exists():
        pytest.skip("not a git checkout")
    tracked = set(subprocess.run(["git", "ls-files", "src/oktopy"], cwd=ROOT, capture_output=True, text=True,
                                 check=True).stdout.splitlines())
    untracked = [f for f in _data_files() if f not in tracked]
    assert not untracked, f"not in git (check .gitignore): {untracked}"


def test_package_data_globs_cover_the_data():
    globs = tomllib.loads((ROOT / "pyproject.toml").read_text())["tool"]["setuptools"]["package-data"]["oktopy"]
    rel = [Path(f).relative_to("src/oktopy").as_posix() for f in _data_files()]
    missed = [f for f in rel if not any(fnmatch.fnmatch(f, g) or fnmatch.fnmatch(f, g.replace("/**/", "/"))
                                        for g in globs)]
    assert not missed, missed
    for f in ("de_cache.R", "de_sensitivity.R", "atlas_template.html", "methods.toml"):
        assert any(fnmatch.fnmatch(f, g) for g in globs), f
