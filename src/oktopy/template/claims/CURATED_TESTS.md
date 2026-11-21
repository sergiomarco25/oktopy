# Curated test settings

Most claims are tested from their contexts (`claims/contexts.tsv`, made by `oktopy contexts` from
`config/contexts.toml`) and their gene signatures (`claims/SIGNATURES.md`). Three small tables, all optional, set the
test of single claims where the rules cannot. Record the reason for every row; re-run `contexts` (overrides) or the
re-tests after a change.

## `claims/contrast_overrides.tsv`: the contrast of one claim

Columns `claim_id`, `contrast`, `note`. Use it when the comparison text maps to the wrong contrast (orientation along a
second variable, a trajectory) or when no contrast in the data matches (`not_available`).

## `claims/state_labels.tsv`: a state's share from the authors' own labels

For claims about the share of an author-defined state (a supertype, subcluster or microglial state) in a dataset whose
nuclei carry that label. The re-test compares the share of the listed labels among the nuclei of `within` per donor x
region (log odds), case vs reference or with the stage. Other datasets keep their usual test (the state's markers, if a
signature exists) unless `other_datasets` is `none`. One row per claim x dataset:

| column | meaning |
|---|---|
| `claim_id` | the claim |
| `dataset` | a dataset whose harmonised obs has the labels |
| `label_column` | the obs column with the labels (kept by the loader, e.g. `supertype_author`) |
| `labels` | the claim's labels, `;`-separated; their nuclei are pooled |
| `within` | the cell classes forming the denominator, `;`-separated (empty: all nuclei of the unit) |
| `measure` | a numeric unit column for a continuous contrast instead of the stage (empty: the stage) |
| `other_datasets` | empty: datasets without the labels keep their usual test (the state's markers); `none`: they cannot test the claim (e.g. a share said unchanged while the state's genes change) |
| `note` | where the labels come from and the rule that chose them |

Choose labels by a rule applied the same way to every claim (e.g. "credible change in >= 80% of the authors'
reference models in both regions"), and name the source table.

## `claims/sex_claims.tsv`: claims about one sex

Claims whose comparison is between the sexes get the contrast `sex`. With `[tests] sex = {column = "sex", contrast =
"<case_vs_control code>"}` in `config/project.toml`, each is tested within the sex it names:

| column | meaning |
|---|---|
| `claim_id` | the claim |
| `sex` | `female` / `male` (the sex the claim says shows the effect; the other sex is reported alongside), or `both` (tested in each sex; reported as the sex that tests the claim hardest: for a no-difference claim the larger change, p doubled for the two tests) |
| `measure` | `composition` (class share, `[tests] sex` contrast, within each sex), `pathology_shift` (mean correlation of all genes with a pathology measure, within the sex; donors permuted), or `not_testable` |
| `column` | for `pathology_shift`: a numeric unit column (empty: the stage) |
| `note` | why this measure tests the claim, or why it cannot be tested |

Without `[tests] sex`, or without a row, a `sex` claim is not testable.
