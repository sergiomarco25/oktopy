# Changelog

## 0.1.0 (unreleased)

First standalone release of the claim-replication engine, until now copied into each project as `src/claimstore/`.

- Engine taken from the reference project's working tree on 2026-11-21 (branch `phase7b-multiverse`, on commit
  2904152, with the uncommitted method-multiverse work: `multiverse.py`, `de_sensitivity.R`,
  `composition_sensitivity.R`, native effect sizes in `retest.py` and `build_store.py`, method registry 0.3.0).
  Rating rules 2.0.0. Docstrings still name the reference project's development phases (7a, 7b) where a choice was
  made.
- Package `oktopy` with one command line, `oktopy`; the R scripts, atlas page, method registry, project template,
  worked example and Claude Code skill ship as package data.
- The engine finds its project from `--project`, `$OKTOPY_PROJECT`, `$OKTOPY_CONFIG` (or the older
  `$CLAIMSTORE_CONFIG`) or the working directory, instead of assuming it sits in `<project>/src/claimstore/`. The
  default environment-variable prefix is `OKTOPY`.
- `oktopy new` (was the reference project's `scripts/63_new_project.py`): projects no longer copy the engine; they
  pin the oktopy version (`requirements.txt`, `[project] oktopy`) and commands warn on a mismatch. New projects get a
  Snakefile with a laptop and a SLURM profile, and inventory columns `n_case` / `n_control`.
- `oktopy example`: two synthetic datasets with planted effects and ten claims of two fictitious papers, end to end
  in minutes. `oktopy run`: every step in order.
- Fixes against the copied engine (outputs of existing projects unchanged unless noted):
  - a claim listing several cell types as a YAML list with no gene or state stopped the store build (`build_store`);
  - the atlas read a field-specific staging column; it now shows the column named by `[claims] staging_field`
    (the page's data gain a `staging` key);
  - the R scripts read counts with `gzip -dc` instead of `zcat` (macOS);
  - the atlas page's always-shown wording is tissue-neutral ("Regions", "Biobank"); a project keeps "brain
    regions" with `[atlas] regions_label`. The brain map and its sections appear only with a project's region
    tables, as before;
  - `multiverse all` and the release's re-rating run the engine as a module, so they work with an installed package.
