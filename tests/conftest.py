"""Test setup: the engine reads a project at import time, so before any test module is imported this makes the worked
example project in a temporary folder and points the engine at it ($OKTOPY_PROJECT), unless one is set already.
The tests then never read or write this repository's own folders. Works with oktopy installed or from the source
tree (src/ is put on the path, for subprocesses too)."""
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
os.environ["PYTHONPATH"] = os.pathsep.join(p for p in (str(SRC), os.environ.get("PYTHONPATH", "")) if p)

EXAMPLE_PAPERS = {"smith2024", "jones2025"}

if not os.environ.get("OKTOPY_PROJECT"):
    from oktopy import example  # noqa: E402  (no engine config read yet)

    _proj = example.make(Path(tempfile.mkdtemp(prefix="oktopy-tests-")) / "example")
    os.environ["OKTOPY_PROJECT"] = str(_proj)


import subprocess  # noqa: E402

import pytest  # noqa: E402


@pytest.fixture(scope="session")
def store_db():
    """The example project's claim store built from its claims alone (no data: every claim Untested); fast, needs
    only pandas, scipy and pyyaml."""
    proj = Path(os.environ["OKTOPY_PROJECT"])
    for step in (["merge"], ["contexts"], ["store"], ["rating"], ["store"]):
        r = subprocess.run([sys.executable, "-m", "oktopy", *step], cwd=proj, capture_output=True, text=True)
        assert r.returncode == 0, (step, r.stdout[-1500:], r.stderr[-1500:])
    from oktopy import config as cfg
    assert cfg.STORE_DB.is_file()
    return cfg.STORE_DB
