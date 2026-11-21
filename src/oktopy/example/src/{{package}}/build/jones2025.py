"""Loader for the example dataset jones2025: synthetic nuclei from 100 donors (50 cases, 50 controls), one region
(front), made in memory by `_synthetic.make` (no download). Its interferon effect is weaker than smith2024's.
"""
from ._synthetic import AUTHOR_TO_CLASS, make  # noqa: F401  (AUTHOR_TO_CLASS: author label -> cell class)


def load():
    return make("jones2025", regions=["front"], n_case=50, n_control=50, seed=2025, ifn_fold=1.6)
