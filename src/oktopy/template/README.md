# {{title}}

{{description}}

A claim store for one field: published conclusions extracted as structured claims (`claims/<paper>.yaml`), re-tested
at donor level in every harmonised dataset, rated for certainty (`claims/RATING.md`) and served to people and AI tools.
Made with `oktopy new` on {{date}}, for oktopy {{oktopy_version}} (`requirements.txt`, `[project] oktopy` in
`config/project.toml`).

## Layout
| Path | What |
|---|---|
| `config/project.toml` | every path, the dataset list, release rules and facts about this field (`oktopy.config`) |
| `config/schema.toml` | the claim schema: required fields and vocabularies; set `region` and `cell_type` for this field |
| `config/contrasts.toml` | case and reference groups of each contrast, and the disease stage |
| `config/harmonise.toml` | the harmonised data schema, cell classes and their marker genes |
| `config/contexts.toml` | how a claim's comparison, cell type and region map to a contrast, cell classes and regions |
| `config/server_instructions.md` | what the MCP server tells AI clients about this store |
| `claims/` | one YAML file of claims per paper (`example2024.yaml` shows the format), findings, reference knowledge |
| `inventory/` | papers screened (`papers.tsv`) and datasets (`datasets.tsv`) |
| `src/{{package}}/` | what is specific to this field: dataset loaders (`build/`), analyses |
| `workflow/` | the steps as a Snakemake workflow, with a laptop and a SLURM profile |
| `tests/` | checks run by CI: the config, the schema, every claim file |
| `.claude/skills/claim-replication/` | the method step by step as a Claude Code skill |

## Commands
Run `oktopy COMMAND` in this folder (`oktopy --help` lists them all), in the order a project runs them:

| Step | Commands |
|---|---|
| claims | `schema`, `merge` (validate and merge claim files), `contexts` (claim -> dataset, regions, cell classes, contrast), `supplements` (list or download the papers' supplementary files), `signatures` (each claim's genes for the re-tests) |
| per dataset | `harmonise DATASET` (loader, schema check, cell classes), `pseudobulk DATASET`, `de-cache DATASET` (edgeR) |
| re-tests | `tests internal DATASET`, `tests cross DATASET`, then `meta` (combine each claim's tests) |
| store | `store`, `rating` (verdict and certainty), `store` again (load the ratings), `serve` / `http` (MCP server), `atlas` (the browsable page, reports/atlas.html; wording in `[atlas]`) |
| robustness | `calibrate`, `power`, `multiverse` (shuffled-label checks and other methods; results/robustness/, results/multiverse/) |
| release | `release VERSION`, `config --show` |

`oktopy run` runs every step from `schema` to `atlas` for the datasets of `[pipeline] datasets`; on a cluster,
`snakemake --profile workflow/profiles/slurm` does the same with one job per dataset and step.
Steps from `contexts` on need pandas; `harmonise` and `pseudobulk` also anndata and scanpy (`oktopy[analysis]`);
`de-cache` needs R with edgeR and data.table (`[envs] r_env` in `config/project.toml`, or Rscript on the PATH).

## First steps
In Claude Code, the `claim-replication` skill walks through every step below and on to a public release, and says
which decisions are yours (papers, datasets, downloads, case and control groups, what stays private).

1. Set `region` and `cell_type` in `config/schema.toml`, then `oktopy schema`; give each region a row in
   `claims/REGIONS.md` (what the code means, the names papers use; checked by `pytest tests`).
2. Add papers to `inventory/papers.tsv` and their claims to `claims/<paper>.yaml`; `pytest tests` validates them.
3. Define contrasts in `config/contrasts.toml`, cell classes and markers in `config/harmonise.toml`, and write one
   loader per dataset in `src/{{package}}/build/` (list them in `[pipeline] datasets`). Large h5ad files:
   `oktopy.h5io.read_csr` / `read_csr_many` stream counts into one matrix without a second copy.
4. Set paths in `config/project.toml` (or `{{env_prefix}}_DATA_ROOT` and friends) and the release rules.

Gene-based claims are tested on the genes each claim names (`signatures`); signatures curated from the papers'
supplementary tables (`supplements --list`, then `supplements`) go in `claims/signatures/<paper>.tsv` and replace
them (`claims/SIGNATURES.md`). A claim whose own data are not harmonised here is still tested in every other dataset.
Claims stay 'Untested' until a dataset re-tests them. The test of single claims can be set where the rules cannot (a
contrast override, a state's share from the authors' own labels, claims about one sex): `claims/CURATED_TESTS.md`.

## Licence
Code: Apache-2.0 (`LICENSE`). Claims and other data and text created here: CC BY 4.0 (`LICENSE-DATA.md`).
