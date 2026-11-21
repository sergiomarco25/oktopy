# Tutorial: the worked example

The example is a complete oktopy project for a made-up field: an "example disease" in an example tissue with two
regions (`front`, `back`) and four cell classes (Epithelial, Immune, Fibroblast, Endothelial). Two fictitious papers,
Smith 2024 and Jones 2025, make ten claims; their two datasets are simulated on the fly (no download), with effects
planted so that every claim has a known answer. Running it shows every step of the method and what it writes.

You need oktopy with the analysis extras and R with edgeR (the conda environment or the Docker image; see the README).

## 1. Make the project

```bash
oktopy example demo
```

```bash
cd demo
```

This is `oktopy new` plus the example's field. Look around:

| file | what to look at |
|---|---|
| `claims/smith2024.yaml`, `claims/jones2025.yaml` | the claims: statement, short quote, comparison, direction, genes, test template |
| `config/schema.toml` | `region` and `cell_type` are the field's own vocabularies |
| `config/contrasts.toml` | `case_vs_control`, a stage contrast (`late_vs_control`) and the stage as a covariate |
| `config/harmonise.toml` | the cell classes and their marker genes |
| `config/contexts.toml` | how "case vs control" or "along disease stage" in a claim becomes a contrast code |
| `src/example_claims/build/smith2024.py` | a dataset loader: returns counts and donor metadata |
| `src/example_claims/build/_synthetic.py` | the simulation, and the table of planted effects |

The project's own checks already pass: `pytest tests` (the claims validate against the schema, the configs agree).

## 2. Run everything

```bash
oktopy run
```

About two minutes on a laptop (most of it edgeR and harmonisation). `oktopy run --list` prints the commands it runs; you
can run any of them yourself, in this order:

**Claims.** `oktopy merge` validates the claim files and merges them into `claims/claims.yaml`; `oktopy contexts`
maps each claim to a dataset, regions, cell classes and a contrast (`claims/contexts.tsv`; before harmonisation every
claim is `no dataset`); `oktopy signatures` collects the genes each claim names (`claims/signatures_all.tsv.gz`).

**Data, per dataset.** `oktopy harmonise smith2024` runs the loader, checks the schema, maps the authors' labels
(`Epi`, `Imm`, ...) to cell classes and checks them against the marker genes, and writes
`data/harmonised/smith2024.h5ad`, a donor table and `smith2024_label_check.tsv` (agreement of each author label with
its marker class). `oktopy pseudobulk smith2024` sums counts per donor x region x cell class
(`data/results/de_cache/smith2024/input/`); `oktopy de-cache smith2024` fits edgeR for every region x class x
contrast and along the stage (`data/results/de_cache/smith2024/de/`, `de_log.tsv`).

**Re-tests.** `oktopy contexts` again (now the datasets' regions are known), then `oktopy tests internal smith2024`
(Smith's claims in Smith's data) and `oktopy tests cross smith2024` (Jones's claims in Smith's data), the same for
Jones, and `oktopy meta` (each claim's tests in the other dataset combined). Results:
`data/results/tests/{internal,cross}/<dataset>.tsv`, one row per claim with its outcome, z, p, donors and the note
for anything not testable.

**Store.** `oktopy store` builds `data/results/store/claim_store.sqlite`, `oktopy rating` rates every claim
(`claims/store/claim_rating.tsv`, with the reason for every certainty step), `oktopy store` loads the ratings, and
`oktopy atlas` writes `reports/atlas.html`.

## 3. Read the results

Open `reports/atlas.html` in a browser, or ask the store directly:

```bash
python -c "import sqlite3; [print(*r) for r in sqlite3.connect('data/results/store/claim_store.sqlite').execute('SELECT claim_id, verdict, certainty FROM claim')]"
```

What you should see, and why:

| claim | says | planted in the data | verdict | certainty |
|---|---|---|---|---|
| smith2024-C01, jones2025-C01 | more immune cells in cases | immune share x2 (epithelial x0.5) | Supported | Moderate, High |
| smith2024-C02 | SPP1 up in immune cells | x3 | Supported | Moderate |
| smith2024-C03, jones2025-C03 | interferon genes up in epithelial cells | x2 (Smith), x1.6 (Jones) | Supported | Moderate, High |
| smith2024-C06 | MKI67 in epithelial cells rises with stage | x1.4 per stage step | Supported | Moderate |
| smith2024-C04 | HIF1A up in endothelial cells | x0.5 (down) | Contradicted | Moderate |
| jones2025-C02 | COL3A1 down in fibroblasts | unchanged | Inconclusive | Low |
| smith2024-C05 | endothelial share does not differ | unchanged | Inconclusive | Low |
| jones2025-C04 | immune cells closer to epithelium (spatial) | no spatial data | Untested | Very low |

Measured on four cores: `oktopy run` takes about 95 seconds (plus Python start-up), the data folder 41 MB.

Things to notice:

- **Both directions of being wrong.** HIF1A (smith2024-C04) is claimed up and is down in the data: Contradicted.
  COL3A1 (jones2025-C02) is claimed down and does not change: not Supported, but Inconclusive rather than No effect,
  because "no effect" needs the confidence interval inside ±0.2 standardised units.
- **No-difference claims need power.** smith2024-C05 says endothelial proportions do not differ. They do not, and
  each dataset's test agrees (outcome reproduced), but equivalence needs a narrow interval: Inconclusive.
- **Untested is not false.** jones2025-C04 is a spatial claim and no dataset here is spatial.
- **Certainty is about the evidence, not the effect.** Two cohorts start at Moderate. A large, precise effect raises
  it a level (Jones's Supported claims: High); a claim loses a level when most of the weight comes from the paper's
  own cohort (dependence), so Smith's claims, weighted towards Smith's larger dataset, stay at Moderate. The two
  no-change claims lose a level for inconsistency or dependence: Low. `claims/store/claim_rating.tsv` gives the
  reason for every step.
- **The donor is the unit.** Every test compares donors; `n_case` and `n_ref` in the test tables count donors, not
  nuclei.

## 4. Change something and see it move

- Add a claim to `claims/jones2025.yaml` (copy one, change the id to `jones2025-C05`, the gene to `MKI67` and the
  comparison to `along disease stage`), then `oktopy run --from merge`.
- Make a planted effect weaker (`ifn_fold` in `src/example_claims/build/jones2025.py`) and rerun from `harmonise`.
- Serve the store to an AI client: `oktopy serve` (MCP over stdio) or `oktopy http`.
- Check the tests themselves: `oktopy calibrate de smith2024` shuffles case and control labels and counts how often
  edgeR still finds something (about 5% is right); `oktopy power` says how much power each test had.

## 5. Your own field

`oktopy new my-field --name my-field --title "..."` makes the same layout with an empty field. The project README
lists the steps; the `claim-replication` skill in `.claude/skills/` walks through them in Claude Code, including
extracting claims from full text, curating gene signatures from supplementary tables and preparing a public release.
