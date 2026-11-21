"""Loader for the example dataset smith2024: synthetic nuclei from 112 donors (56 cases, 56 controls), two regions
(front, back), made in memory by `_synthetic.make` (no download). A real loader reads the authors' files here and
selects donors (state the rule in this docstring); see the oktopy tutorial.
"""
from ._synthetic import AUTHOR_TO_CLASS, make  # noqa: F401  (AUTHOR_TO_CLASS: author label -> cell class)


def load():
    return make("smith2024", regions=["front", "back"], n_case=56, n_control=56, seed=2024, ifn_fold=2.0)
