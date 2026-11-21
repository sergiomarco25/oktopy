# Config reference

A project is a folder with `config/project.toml`. oktopy finds it from `--project DIR`, else `$OKTOPY_PROJECT`,
else the nearest folder at or above the working directory that holds `config/project.toml`
(`oktopy config --show` prints what it found). Every file below is made by `oktopy new` with comments; this page
says what each setting does and which command reads it.

| file | what | read by |
|---|---|---|
| `config/project.toml` | paths, pipeline, release rules, facts about the field | every command |
| `config/schema.toml` | the claim schema: required fields, vocabularies, field descriptions | `schema`, `merge`, project tests |
| `config/contrasts.toml` | which donors are case and reference in each contrast; the stage | `pseudobulk`, `de-cache`, `tests` |
| `config/harmonise.toml` | the harmonised data schema, cell classes, marker genes, labelling rules | `harmonise`, `pseudobulk`, `tests` |
| `config/contexts.toml` | how claim text maps to contrasts, cell classes and regions | `contexts`, `store` |
| `config/server_instructions.md` | what the MCP server tells AI clients | `serve`, `http` |
| `config/methods.toml` (optional) | methods of the registry this project switches off | `registry`, `multiverse` |

## `config/project.toml`

### `[project]`
| key | meaning |
|---|---|
| `name`, `title`, `description`, `keywords` | how releases, the atlas, the server and linked-data exports name the store |
| `server_name` | the MCP server's name |
| `licence` | licence statement of releases |
| `env_prefix` | prefix of the environment variables that override paths (`<PREFIX>_DATA_ROOT`, `<PREFIX>_STORE_DB`, ...) |
| `oktopy` | the oktopy version the project was made and checked with; commands warn when another version runs |

### `[paths]` and `[envs]`
| key | meaning |
|---|---|
| `data_root` | large files, out of git (default `data`); relative paths are relative to the project folder |
| `work` | where `raw/`, `harmonised/`, `results/` and `scratch/` go (default `{data_root}`) |
| `python` | another Python for steps run by a separate interpreter (the release's re-rating) |
| `r_env` | an R environment folder whose `bin/Rscript` has edgeR and data.table (empty: `Rscript` on the PATH) |

Values may use `{repo}` (the project folder) and `{data_root}`, `~` and environment variables. Each can be overridden
by `<PREFIX>_<NAME>` (for example `MY_FIELD_DATA_ROOT=/fast/disk/my-field`). On a cluster, put `data_root` on a shared
file system all nodes can read; keep the project folder (config, claims, code) in git.

### `[pipeline]`
| key | meaning |
|---|---|
| `loaders` | the package holding one loader module per dataset (default `<package>.build`) |
| `datasets` | the datasets `oktopy run` and the Snakefile build and test |
| `default_mem_gb`, `mem_gb` | memory of per-dataset cluster jobs (Snakefile), default and per dataset |

### `[claims]`
| key | meaning |
|---|---|
| `paper_to_dataset` | papers whose claims were made on another paper's cohort: paper id -> dataset id |
| `sample_type` | how the store describes samples ("post-mortem tissue", "biopsy", "organoid") |
| `heritability_label` | subject of genetic-enrichment claims |
| `stage_available_column`, `staging_field` | the store's dataset columns for "any donor staged" and the staging description (cohort file) |

### `[tests]`
| key | meaning |
|---|---|
| `gene_aliases` | gene symbols renamed since publication: `old = "new"` (MARCH*/SEPT* renames are built in) |
| `concordance` | contrast codes comparing two contrasts' fold changes: `code = {contrasts = ["a", "b"], template = "name"}` |
| `genetic_test` | `"module:function"` of the project's package testing `genetic_enrichment` claims (empty: not testable) |
| `signatures_per_class` | apply a curated signature per cell class when its rows name classes (recommended: true) |
| `min_score_genes` | smallest gene set tested with the donor score (default 2) |
| `unit_columns` | cell metadata the composition test reads (default: the harmonised schema's donor columns) |
| `sex` | `{column = "sex", contrast = "<code>"}`: claims with the contrast `sex` are tested within each sex (`claims/sex_claims.tsv`) |

### `[rating]`
| key | meaning |
|---|---|
| `effect_templates`, `pattern_templates` | test templates of each rating family (`claims/RATING.md` section 2) |
| `enriched_sort` | datasets of sorted nuclei or cells, whose composition is not tissue composition |
| `shared_banks`, `shared_bank_caveat` | datasets that may share donors (a caveat on the rating, not merged) |
| `not_assessed` | what the rating does not assess, stated on every rating |

### `[atlas]` (all optional)
`title`, `headline`, `lead`, `crumb` (page wording; default the project's title and description), `disease`,
`case_label` (names for the case group), `scope` (what the records are about, for "Ask the atlas"), `examples`,
`ask_placeholder`, `stop_words`, `region_groups` / `region_group_of` (columns of the data landscape), `key_findings`
(`[{label, title, text, href}]` cards), `graph_start` (node ids), `template` (a page of the project's own), `extras`
(`"module:function"` adding tables), `out` (default `reports/atlas.html`).

### `[server]`
`instructions` (a file next to `project.toml`) and `examples` (ids quoted in the tool descriptions: a claim, gene,
cell type, region, dataset, paper of this store).

### `[release]`
| key | meaning |
|---|---|
| `public_repository` | where releases and the public atlas are published |
| `private_datasets`, `private_papers` | kept out of every public export (unpublished data) |
| `separable_findings` | findings that keep their public datasets when private ones are removed |
| `rederived_findings` | findings recomputed without private data: `id = {table, script}` |
| `private_modules` | names (of programmes, states) that must not appear in public exports |

## `config/schema.toml`
`required` and `list_fields` (fields that may hold a list), `[vocab]` (allowed values; set `region` and
`cell_type` for the field and keep the others), `[text]` and `[finding_text]` (field descriptions, written into the
JSON Schemas). After a change: `oktopy schema`, and give every region code a row in `claims/REGIONS.md`.

## `config/contrasts.toml`
`[contrasts.<code>]` with `case` and `ref` tables of conditions on the unit metadata; a unit belongs to a side when
all conditions hold: `column = "value"`, `column = ["a", "b"]`, `column = {contains = "x"}`, or a stage range
`stage = {min = 2, max = 4}`. `[stage]`: `name`, `columns` (the first with values is the stage), `group_column` and
`zero_group` (that group counts as stage 0), `contrast` (the code claims use for a continuous association) and
`label`. Top-level lists are named groups.

## `config/harmonise.toml`
`obs` (cell metadata every harmonised dataset carries), `donor_extras` (more donor-level columns for the unit
table), `regions`, `diseases`, `cell_classes` (with `Unknown`), `sorted_classes` (classes tested on sorted nuclei in
datasets mixing sorted and unsorted nuclei), `[markers]` (canonical genes per class, at least 2 detected), and
`[labelling]`: a two-step call (lineage, then neuron subtype: `neuron_subtypes`, `pan_neuronal`, `subtype_min_z`,
`precedence`; leave empty for a one-step call), `min_raw`, `[labelling.region_only]` (classes plausible only in some
regions), `[labelling.require_detect]` (classes a cluster takes only when a gene is detected in enough nuclei).

## `config/contexts.toml`
`null_pattern` (regex recognising no-difference statements), `rules` (ordered `[regex, contrast]` pairs on the
lower-cased comparison text), `[cell_class_expand]` (claim cell type -> tested classes), `[region_map]` (claim region
-> dataset regions), `[contrast_labels]` (code -> `[kind, case group, reference group]`; kinds: case_control, stage,
continuous, covariate, within, not_available, unassigned).

## Files that set tests of single claims
`claims/contrast_overrides.tsv` (claim_id, contrast, note), `claims/state_labels.tsv` (a state's share from the
authors' labels) and `claims/sex_claims.tsv` (claims about one sex): `claims/CURATED_TESTS.md` in every project.
Curated gene signatures: `claims/signatures/<paper>.tsv` (`claims/SIGNATURES.md`).

## Environment variables
| variable | meaning |
|---|---|
| `OKTOPY_PROJECT` | the project folder (`oktopy --project DIR` sets it) |
| `OKTOPY_CONFIG` | a config file other than `<project>/config/project.toml` |
| `<PREFIX>_DATA_ROOT`, `<PREFIX>_WORK`, `<PREFIX>_PYTHON`, `<PREFIX>_R_ENV` | override `[paths]` and `[envs]` |
| `<PREFIX>_STORE_DB` | the store the server opens |
| `<PREFIX>_MCP_READONLY=1`, `<PREFIX>_MCP_ALLOW_WRITE=1` | stdio server without `add_question`; HTTP server with it |
| `<PREFIX>_HTTP_HOST`, `PORT`, `<PREFIX>_PUBLIC_URL`, `<PREFIX>_MCP_TOKEN`, `<PREFIX>_ALLOWED_ORIGINS`, `<PREFIX>_RATE_PER_MIN`, `<PREFIX>_SQL_SECONDS` | the HTTP server (`oktopy.mcp_http`) |
