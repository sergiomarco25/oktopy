"""oktopy: re-test published claims of a field at donor level in every harmonised dataset, rate them and serve them.

Usage: oktopy [--project DIR] COMMAND [ARGS...]      (oktopy COMMAND --help for a command's options)

Start a project
  new         DIR --name NAME --title TITLE: a new project for your field (its checks run at once)
  example     DIR: a worked example project with two small synthetic datasets and two papers' claims
  run         every step below in order for the project's datasets ([pipeline] datasets); --from / --to STEP
Claims
  schema      write claims/schema/*.json from config/schema.toml
  merge       validate the per-paper claim files and merge them (claims/claims.yaml, claims.tsv, validation.tsv)
  contexts    map every claim to its dataset, regions, cell classes and contrast (claims/contexts.tsv)
  supplements [--list] [PAPER...]: list (--list) or download the papers' supplementary files (raw/supplements/)
  signatures  claim gene signatures for the re-tests: curated claims/signatures/*.tsv plus the genes each claim names
Data, per dataset
  harmonise   DATASET: run the dataset's loader, check the schema, label cell classes, write harmonised files
  relabel     DATASET...: re-run labelling checks, donor table and pseudobulk after a change to the rules
  pseudobulk  DATASET: export counts per region x cell class and the unit table for edgeR
  de-cache    DATASET...: edgeR fits per region x class x contrast and for the stage (R: [envs] r_env)
  de-log      DATASET...: rebuild de_log.tsv from written fit files
Re-tests, rating, store
  tests       internal|cross DATASET [gwas|nogwas]: re-test claims in one dataset (results/tests/)
  meta        combine every claim's tests across datasets (results/meta/claim_meta.tsv)
  store       build the claim store (results/store/claim_store.sqlite, claims/store/)
  rating      rate every claim of the store (claims/store/claim_rating.tsv); --db / --out
  atlas       the browsable atlas page of the store (reports/atlas.html); --public --db STORE for a release
Robustness
  registry    validate the method registry (oktopy/methods.toml); --md OUT writes it as a page
  calibrate   de|sets|scores|composition|state|summary DATASET...: null checks of the current tests (results/robustness/)
  power       detectability per claim x dataset and reasons for Inconclusive, from the store (results/robustness/)
  multiverse  sets|fit|null|tests|summary: claims re-tested with the registry's sensitivity methods (results/multiverse/)
Release and serving
  release     build a public release (VERSION ...);  diff OLD NEW: per-claim changes;  linked-data RELEASE_DIR
  audit       PAGE: check a public atlas page for private names
  serve       MCP server over stdio;  http: MCP + REST server over HTTP
  config      print paths (--show, --env, NAME)

The project folder is --project DIR, else $OKTOPY_PROJECT, else the nearest folder at or above the working directory
holding config/project.toml. Steps from `contexts` on need pandas; harmonise, relabel and pseudobulk also anndata and
scanpy (pip install "oktopy[analysis]"); de-cache needs R with edgeR and data.table.
"""
import importlib
import os
import sys

from oktopy import __version__

COMMANDS = {"schema": "schema", "merge": "claims", "contexts": "contexts", "signatures": "signatures",
            "supplements": "supplements", "harmonise": "harmonise", "relabel": "harmonise:relabel_main",
            "pseudobulk": "pseudobulk", "de-cache": "de_cache", "de-log": "de_cache:rebuild_main", "tests": "retest",
            "meta": "claim_meta", "store": "build_store", "rating": "rating", "registry": "registry",
            "calibrate": "calibration", "power": "power", "multiverse": "multiverse", "atlas": "atlas",
            "release": "release", "serve": "mcp_server", "http": "mcp_http", "config": "config",
            "diff": "release_diff", "linked-data": "linked_data", "audit": "audit",
            "new": "newproject", "example": "example", "run": "runner"}
# commands that make or drive a project rather than read one: no version check against the project's pin
NO_PIN_CHECK = {"new", "example", "config"}


def _pin_warning():
    """Warn when the project pins another oktopy version ([project] oktopy in config/project.toml)."""
    from oktopy import config as cfg
    pin = str(cfg.PROJECT.get("oktopy", "")).strip()
    if pin and pin != __version__:
        print(f"oktopy: this project was made for oktopy {pin}, installed is {__version__} "
              f"(pin it: pip install oktopy=={pin}, or update [project] oktopy after checking the changes)",
              file=sys.stderr)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--version"]:
        print(f"oktopy {__version__}")
        return 0
    if argv[:1] == ["--project"]:
        if len(argv) < 2:
            print("--project needs a folder", file=sys.stderr)
            return 2
        os.environ["OKTOPY_PROJECT"] = os.path.abspath(argv[1])
        argv = argv[2:]
    elif argv and argv[0].startswith("--project="):
        os.environ["OKTOPY_PROJECT"] = os.path.abspath(argv[0].split("=", 1)[1])
        argv = argv[1:]
    if not argv or argv[0] in ("-h", "--help") or argv[0] not in COMMANDS:
        print(__doc__)
        if argv and argv[0] not in ("-h", "--help"):
            print(f"unknown command {argv[0]!r}", file=sys.stderr)
        return 0 if argv and argv[0] in ("-h", "--help") else 2
    if argv[0] not in NO_PIN_CHECK:
        _pin_warning()
    name, _, fn = COMMANDS[argv[0]].partition(":")
    mod = importlib.import_module(f"oktopy.{name}")
    sys.argv = [f"oktopy {argv[0]}", *argv[1:]]
    out = getattr(mod, fn or "main")()
    return out if isinstance(out, int) else 0


if __name__ == "__main__":
    sys.exit(main())
