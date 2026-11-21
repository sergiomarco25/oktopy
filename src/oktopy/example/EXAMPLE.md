# The oktopy worked example

Two fictitious papers (Smith 2024, Jones 2025) make ten claims about an example disease; their datasets are
simulated by `src/{{package}}/build/_synthetic.py` with effects planted in the case donors, so every claim has a
known answer. Run everything with `oktopy run` (about two minutes), then open `reports/atlas.html`.

| claim | says | planted in the data | expected verdict |
|---|---|---|---|
| smith2024-C01, jones2025-C01 | more immune cells in cases | immune share x2 (epithelial x0.5) | Supported |
| smith2024-C02 | SPP1 up in immune cells | x3 | Supported |
| smith2024-C03, jones2025-C03 | interferon genes up in epithelial cells | x2 (Smith), x1.6 (Jones) | Supported |
| smith2024-C06 | MKI67 in epithelial cells rises with stage | x1.4 per stage step | Supported |
| smith2024-C04 | HIF1A up in endothelial cells | x0.5 (down) | Contradicted |
| jones2025-C02 | COL3A1 down in fibroblasts | unchanged | Inconclusive |
| smith2024-C05 | endothelial share does not differ | unchanged | Inconclusive (equivalence needs more donors) |
| jones2025-C04 | immune cells closer to epithelium (spatial) | no spatial data | Untested |

Certainty and the reason for each step: `claims/store/claim_rating.tsv`. The walk-through: docs/tutorial.md in
the oktopy repository.
