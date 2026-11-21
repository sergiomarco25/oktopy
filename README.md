# oktopy

**Which published single-cell findings hold up?** oktopy turns the conclusions of a field's papers into structured
claims, re-tests every claim at donor level in its own dataset and in every other harmonised dataset, pools the tests,
rates each claim (a verdict and a GRADE-like certainty) and serves the result as a SQLite claim store, a knowledge
graph, a browsable atlas page and an MCP server for AI tools.

It is the engine behind two public claim stores, for Parkinson's disease
([pd-sc-claims](https://github.com/sergiomarco25/pd-sc-claims)) and Alzheimer's disease
([ad-sc-claims](https://github.com/sergiomarco25/ad-sc-claims)). Nothing in the engine names a disease: a field is
config files, one loader per dataset and the claims.

```
papers ──> claims/<paper>.yaml ──┐
                                 ├─> contexts ─> re-tests in every dataset ─> meta-analysis ─> rating ─> store, atlas, MCP
datasets ─> loaders ─> harmonise ─> pseudobulk ─> edgeR DE cache ──┘
```

## Install

oktopy needs Python 3.11+ and, for the differential-expression step, R with edgeR and data.table.

```bash
git clone https://github.com/sergiomarco25/oktopy && cd oktopy
```

Then one of:

```bash
mamba env create -f envs/oktopy.yml && mamba activate oktopy && pip install -e .
```

```bash
docker build -t oktopy .
```

```bash
pip install -e ".[analysis,test]"
```

The conda environment and the Docker image have everything (Python, R with edgeR, limma, DESeq2, Snakemake). The
plain `pip` route gives the Python side only: claims, store, rating, atlas and server work, and the data steps need R
with edgeR and data.table on the PATH (or `[envs] r_env` in the project config). Check the install with
`oktopy --version`.

## Quickstart: the worked example (about 2 minutes on a laptop)

```bash
oktopy example demo
```

```bash
cd demo && oktopy run
```

`oktopy example` makes a project with two synthetic single-nucleus datasets (no download) and ten claims from two
fictitious papers whose true answers are planted in the data. `oktopy run` harmonises both datasets, builds pseudobulk
and the edgeR DE cache, re-tests every claim in both datasets, pools the tests, rates every claim and writes the store
and the atlas. Open `demo/reports/atlas.html` in a browser. The [tutorial](docs/tutorial.md) walks through every step
and what it wrote, and lists the verdicts you should see (Supported, Contradicted, Inconclusive, Untested).

With Docker:

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/work oktopy example demo
```

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD/demo":/work oktopy run
```

## Start a project for your field

```bash
oktopy new my-field --name my-field --title "My field single-cell claim store"
```

The new project's own checks run at once (`pytest tests` in it). Then, in order: the field's regions and cell types
(`config/schema.toml`), papers and their claims (`claims/<paper>.yaml`), contrasts (`config/contrasts.toml`), cell
classes and markers (`config/harmonise.toml`), one loader per dataset (`src/<package>/build/<dataset>.py`), and
`oktopy run`. The project's README lists the steps; the [config reference](docs/config.md) explains every setting.

In [Claude Code](https://claude.com/claude-code), the `claim-replication` skill that every project gets
(`.claude/skills/claim-replication/`) walks through the whole method, including claim extraction from full text, and
says which decisions are yours.

## Where it runs

| Where | How |
|---|---|
| Laptop or workstation | `oktopy run` in the project folder (one step after another), or `snakemake` (4 cores by default) |
| SLURM cluster | `snakemake --profile workflow/profiles/slurm` in the project folder: one job per dataset and step; set partition and QOS in the profile, memory per dataset in `[pipeline] mem_gb` |
| Container, cloud VM | the Docker image (`oktopy` is its entry point); on a cluster, `apptainer build oktopy.sif apptainer.def` |
| Serving | `oktopy serve` (MCP over stdio) or `oktopy http` (MCP and REST over HTTP) on a built store |

[docs/running.md](docs/running.md) has the details: profiles, paths on shared file systems, memory, containers.

## Commands

`oktopy --help` lists them all; `oktopy --project DIR COMMAND` runs one in another folder.

| Step | Commands |
|---|---|
| start | `new`, `example`, `run` (every step below in order) |
| claims | `schema`, `merge`, `contexts`, `supplements`, `signatures` |
| per dataset | `harmonise DATASET`, `pseudobulk DATASET`, `de-cache DATASET` |
| re-tests | `tests internal DATASET`, `tests cross DATASET`, `meta` |
| store | `store`, `rating`, `atlas`, `serve`, `http` |
| robustness | `registry`, `calibrate`, `power`, `multiverse` |
| release | `release VERSION`, `diff`, `linked-data`, `audit`, `config` |

## Documentation

- [Tutorial](docs/tutorial.md): the worked example step by step.
- [Concepts](docs/concepts.md): claims and their schema, test templates, outcomes, pooling, the rating, robustness
  and the method multiverse, the store and how it is served.
- [Config reference](docs/config.md): every file and table of a project's `config/`.
- [Running](docs/running.md): laptop, SLURM, containers, cloud.
- [Rating rules](src/oktopy/template/claims/RATING.md): the verdict and certainty rules, versioned.
- [The skill](src/oktopy/skill/claim-replication/SKILL.md): the method as a Claude Code skill.

## Develop

```bash
pip install -e ".[all]" && pytest
```

The tests build the worked example in a temporary folder; the end-to-end test runs it fully when R with edgeR is
available (CI requires it). The engine must name no field (`tests/test_engine_boundary.py`); a change that could
move a verdict needs a test and a note in `CHANGELOG.md`, and projects pin the oktopy version they were checked with
(`requirements.txt`, `[project] oktopy`).

## Licence and citation

Code: Apache-2.0 (`LICENSE`). Projects made with `oktopy new` license their claims, store and text under CC BY 4.0
(their `LICENSE-DATA.md`). Cite oktopy with `CITATION.cff`.
