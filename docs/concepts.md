# Concepts

What oktopy does with a field's papers and datasets, and what each result means. The [tutorial](tutorial.md) shows
the same steps on the worked example; the [config reference](config.md) lists every setting.

## 1. Claims

A **claim** is one conclusion a paper draws from its own data, written so that another dataset can test it: a group
difference, an association, a subpopulation, a gene's change in a cell type, a pathway's change, a relation to
pathology or a clinical measure. Each paper's claims are a YAML list in `claims/<paper_id>.yaml`.

| field | meaning |
|---|---|
| `id`, `paper` | `<paper_id>-C<nn>`, numbered in order of appearance; the paper id is the file name |
| `statement` | the conclusion in one self-contained sentence (cell type, region, groups, direction) |
| `quote`, `location` | the authors' words supporting it (at most 15 words) and where (figure, section) |
| `centrality` | headline (abstract or title), main (Results heading or main figure) or secondary |
| `modality`, `region`, `cell_type`, `cell_state` | what the claim is about, from the project's vocabularies; several values `;`-separated (`Immune; Epithelial`) |
| `claim_type` | composition, de_gene, de_pathway, state_existence, pathology_association, ... |
| `direction` | up, down, enriched, depleted, exists, associated_pos, associated_neg, or NA |
| `comparison` | the groups or variable in the authors' words ("case vs control", "along disease stage") |
| `genes`, `gene_source` | the genes the claim names, and where its full gene list is |
| `reported_stat`, `unit_of_analysis`, `method` | what the authors reported and how they tested it (cells, donors, ...) |
| `testable_with`, `test_template` | the data a re-test needs and the template that runs it (section 3) |

The project's `config/schema.toml` holds the required fields and the vocabularies; `oktopy schema` writes the JSON
Schemas that `pytest tests` and `oktopy merge` check every file against. Only `region` and `cell_type` are the
field's own; the other vocabularies are shared by every oktopy project, so claims stay comparable across fields.
A claim that asserts *no* difference (direction NA, "do not differ", "unchanged") is a claim too: `null_pattern` in
`config/contexts.toml` recognises it and the rating judges it on equivalence.

**Findings** (`claims/findings.yaml`) are claims made by the project's own analyses, rated the same way;
**reference knowledge** (`claims/background.yaml`) records evidence of other kinds (genetics, imaging, animal work)
that can support a claim; **questions** (`claims/questions.yaml`) are what users asked of the store.

## 2. Contexts: from a claim's words to a test

`oktopy contexts` maps every claim to the data that can test it (`claims/contexts.tsv`):

- **contrast**: the comparison text is matched against ordered rules in `config/contexts.toml` (`rules`); a contrast
  code is either a group contrast of `config/contrasts.toml` (which donors are case and which are reference), the
  stage (a continuous covariate), `within` (a structure inside one group) or `not_available` (no dataset can make
  that comparison). `claims/contrast_overrides.tsv` sets it per claim when the rules cannot.
- **cell classes and regions**: claim cell types expand to the harmonised cell classes (`[cell_class_expand]`;
  `All` is the pseudobulk of all cells, `ALL`), claim regions to dataset regions (`[region_map]`).
- **status**: `ok`, `no dataset` (the claim's own data are not harmonised here: it is still tested in every other
  dataset), `not_testable` (template `not_testable`), `not_available`, or `no contrast` (no rule matched: fix the
  rules or add an override; such a claim is never tested).

Direction is always relative to the contrast's case group. Check the orientation of every claim whose comparison
mentions a second variable ("along decline" and "along function" have opposite signs).

## 3. Data: harmonise, pseudobulk, differential expression

- **Harmonise** (`oktopy harmonise DATASET`): the dataset's loader (`src/<package>/build/<dataset>.py`) returns raw
  integer counts and the schema's cell metadata (donor, sample, region, disease, sex, age, pmi, stage, the authors'
  cell label). Harmonisation maps the authors' labels to the project's cell classes and checks them against marker
  genes (`harmonised/<dataset>_label_check.tsv`: agreement, and each label's most common marker class), or labels
  clusters by markers when there are no author labels. It stops on missing donors or regions, unknown regions,
  diagnoses or cell classes, and non-integer counts.
- **Pseudobulk** (`oktopy pseudobulk DATASET`): counts summed per donor x region x cell class (and `ALL`), units
  with at least 10 nuclei, and a unit table with one column per contrast (case, ref or empty).
- **DE cache** (`oktopy de-cache DATASET`, R): edgeR quasi-likelihood fits for every region x class x contrast with
  at least 3 case and 3 reference units, and along the stage; covariates age, sex and PMI when known for every unit
  and estimable. Each fit is a table of genes (`results/de_cache/<dataset>/de/`), logged in `de_log.tsv`.

**The donor is the unit.** Every test compares donors (or donor x region units), never cells; cell-level p values
treat thousands of nuclei from a few people as independent and are not used.

## 4. Re-tests and outcomes

`oktopy tests internal DATASET` tests the claims whose source data the dataset is (internal replication);
`oktopy tests cross DATASET` tests every other claim in it, with the claim's regions retargeted to the dataset's
(cross-dataset replication). The test template decides how:

| template | test | unit |
|---|---|---|
| `composition_donor` | centred log-ratio share of the class per unit; Mann-Whitney case vs reference, or Spearman with the stage | donor x region |
| `de_pseudobulk` (one gene) | the gene's edgeR result, FDR < 0.05 | donor x region (edgeR fit) |
| `de_pseudobulk`, `geneset_shift`, `state_signature` (2+ genes) | donor signature score: mean of the genes' z-scored log CPM, each with its expected sign; Mann-Whitney or Spearman | donor x region |
| `pathology_correlation` | the same, along the stage or a measure | donor x region |
| `state_signature` with `within` | state existence: the genes' co-expression in nuclei against random genes, with donor spread | nuclei, checked across donors |
| `genetic_enrichment` | the field's own test (`[tests] genetic_test`, a function of the project's package) | |
| labels, sex (curated) | a state's share from the authors' own labels (`claims/state_labels.tsv`); claims about one sex (`claims/sex_claims.tsv`) | donor x region |

Regions and cell classes a claim covers are tested separately and combined by a weighted Stouffer z (weights the
square root of the units). Every result row has an **outcome**:

| outcome | meaning |
|---|---|
| `reproduced` | significant in the claimed direction (FDR < 0.05 for a single gene, p < 0.05 for a set or share) |
| `direction_only` | the claimed direction, not significant |
| `not_reproduced` | the opposite direction, not significant (or no expected direction and not significant) |
| `contradicted` | significant in the opposite direction |
| `not_testable` | the data, groups or genes are missing (the reason is in `note`) |

A no-difference claim is `reproduced` when its test is not significant and `contradicted` when it is; the rating then
asks for equivalence, which needs far more donors than detecting a difference.

**Gene signatures.** Gene-based claims are tested on the genes they name (`oktopy signatures`). A paper's full
gene list, curated from its supplementary tables into `claims/signatures/<paper>.tsv` (rules in
`claims/SIGNATURES.md`), replaces the named genes. Curate before trusting gene-level verdicts: three to five named
genes can disagree with the paper's full table.

## 5. Pooling across datasets

`oktopy meta` combines each claim's tests in every other dataset (`results/meta/claim_meta.tsv`): a weighted
Stouffer z and I² heterogeneity, by scope (all other datasets, same region, other regions, datasets published after
the paper) and a robustness class (robust, supported, contradicted, inconsistent, not_replicated, insufficient).

## 6. Rating: verdict and certainty

`oktopy rating` rates every claim from its informative tests (outcome other than `not_testable`), GRADE-style. The
rules are versioned in `claims/RATING.md`; in short:

- **Effect family** (composition, DE genes, gene-set and state shifts): each test's oriented z becomes an approximate
  standardised mean difference, pooled across datasets by random effects. **Pattern family** (state existence,
  genetic enrichment, concordance): detection across datasets and a combined z.
- **Verdict**: Supported (CI excludes 0 on the claimed side), Contradicted (on the opposite side), No effect (CI inside
  ±0.2), Inconclusive (anything else), Untested (no informative test: not the same as false). A no-difference claim
  is "Supported (no effect)" when the CI lies inside ±0.2.
- **Certainty**: starts from the number of independent cohorts (3+ High, 2 Moderate, 1 Low), then loses a level for
  each serious problem: risk of bias (sorted nuclei), inconsistency (I², disagreeing datasets), indirectness (other
  region, partial proxy, non-human evidence), imprecision (wide CI, fewer than 200 donor units; two levels below 30
  or a single test), dependence (most weight from the paper's own cohort). Large effects and support from other
  evidence types can raise it.

Every rating carries its reasons, the regions tested and those where the evidence held, and caveats (datasets that
may share donors).

## 7. Robustness and the method multiverse

The verdict rests on one primary method per task. oktopy checks that choice in three ways, none of which changes a
verdict today:

- **Method registry** (`oktopy registry`, `src/oktopy/methods.toml`): for every task (DE case-control, DE along the
  stage, small and large gene sets, composition, author labels, ...) the methods that could run, ranked into tiers
  (primary, sensitivity, exploratory, excluded) by recommended practice (the single-cell best-practices book),
  donor-level units, calibration and benchmarks, each with its evidence and reason. A project can switch methods off
  in its own `config/methods.toml`, never promote one.
- **Calibration** (`oktopy calibrate`): each test on shuffled case/reference labels, where no difference is real; a
  calibrated test gives p < 0.05 in 2.5-7.5% of null tests (`results/robustness/calibration_summary.tsv`).
- **Power** (`oktopy power`): for every informative test, the power it had to detect a trivial (0.2), important (0.5)
  and the pooled effect, and why each Inconclusive claim is inconclusive (underpowered, heterogeneous, small if any).
- **Multiverse** (`oktopy multiverse`): every claim re-tested with the registry's sensitivity methods in place of its
  primary (limma-voom, DESeq2, Wilcoxon or Spearman per gene; camera and the covariate-adjusted donor score for sets;
  propeller for composition). A method counts in a dataset only after passing the same shuffled-label check there.
  `results/multiverse/agreement.tsv` lists, per claim and dataset, whether the methods agree with the primary.

## 8. The store and how it is served

`oktopy store` builds `results/store/claim_store.sqlite`: claims with their structure (subject, predicate, contrast,
context) and rating, datasets and donor groups, every test (`claim_dataset`), signatures, papers, relations between
claims, consensus statements, reference knowledge, questions and vocabularies. Run it, then `oktopy rating`, then
`oktopy store` again to load the ratings.

- **Knowledge graph** (`oktopy.kg`): the store as nodes (claims, genes, cell types, regions, datasets, papers) and
  edges, exported with a release (GraphML, JSON) and traversed by the server.
- **Atlas** (`oktopy atlas`): one self-contained HTML page (`reports/atlas.html`, works offline) to browse claims,
  ratings, tests, datasets, genes and the graph, with an "Ask the atlas" box. Wording comes from `[atlas]`.
- **MCP server**: `oktopy serve` (stdio) and `oktopy http` (MCP Streamable HTTP plus REST and OpenAPI) give AI tools
  read-only access: `overview`, `search`, `search_claims`, `get_claim`, `get_dataset`, `get_statement`, `gene`,
  `neighbours`, `connect`, `sql`, `schema`, `list_questions` (and `add_question` when writing is allowed).

## 9. Releases and private data

`oktopy release VERSION` writes a public release (`results/releases/v<VERSION>/`): the store without private
datasets and papers (`[release]` in `config/project.toml`), every claim re-rated on what remains, TSV tables, the
knowledge graph, schema.org ClaimReview JSON-LD, an RO-Crate description, the rating rules, a changelog against the
previous release and checksums. `oktopy atlas --public --db <release store>` makes the public atlas and
`oktopy audit PAGE` checks it for private names. Unpublished data stay listed under `[release]` and out of every
export until they are published.
