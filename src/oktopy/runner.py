"""Every step of a project in order, for its datasets ([pipeline] datasets in config/project.toml), one process per
step. For a laptop or a single machine; on a cluster use the project's workflow/Snakefile with the slurm profile.

Usage: oktopy run [--datasets A,B] [--from STEP] [--to STEP] [--list] [--dry-run]

Steps: schema, merge, contexts, signatures, harmonise, pseudobulk, de-cache, contexts-data (contexts again, now that
the datasets' regions are known), tests (internal and cross in every dataset), meta, store, rating, store-rated
(the store again, with the ratings), atlas. Per-dataset steps run for every dataset before the next step starts.
A failing step stops the run with its exit code; --from STEP resumes there.
"""
import argparse
import shutil
import subprocess
import sys
import time

from oktopy import config as cfg

STEPS = ["schema", "merge", "contexts", "signatures", "harmonise", "pseudobulk", "de-cache", "contexts-data", "tests",
         "meta", "store", "rating", "store-rated", "atlas"]
PER_DATASET = {"harmonise", "pseudobulk", "de-cache"}


def commands(step, datasets):
    """The oktopy commands (argument lists) of one step."""
    if step in PER_DATASET:
        return [[step, ds] for ds in datasets]
    if step == "contexts-data":
        return [["contexts"]]
    if step == "tests":
        return [["tests", mode, ds] for ds in datasets for mode in ("internal", "cross")]
    if step == "store-rated":
        return [["store"]]
    return [[step]]


def check_r():
    from oktopy.de_cache import rscript
    r = rscript()
    if shutil.which(r) is None and not cfg.R_ENV.joinpath("bin", "Rscript").is_file():
        return (f"no Rscript found (looked for {r}): install R with edgeR and data.table (envs/oktopy.yml or the "
                "Docker image), or set [envs] r_env in config/project.toml to an R environment folder")
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(prog="oktopy run", description=__doc__.split("\n\n")[0])
    ap.add_argument("--datasets", help="comma-separated datasets (default: [pipeline] datasets)")
    ap.add_argument("--from", dest="start", choices=STEPS, default=STEPS[0])
    ap.add_argument("--to", dest="stop", choices=STEPS, default=STEPS[-1])
    ap.add_argument("--list", action="store_true", help="print the commands and stop")
    ap.add_argument("--dry-run", action="store_true", help="same as --list")
    a = ap.parse_args(sys.argv[1:] if argv is None else argv)
    datasets = [d for d in (a.datasets.split(",") if a.datasets else cfg.DATASETS) if d]
    steps = STEPS[STEPS.index(a.start): STEPS.index(a.stop) + 1]
    if not steps:
        ap.error("--from comes after --to")
    if not datasets and any(s in PER_DATASET or s == "tests" for s in steps):
        print("no datasets: list them in [pipeline] datasets of config/project.toml (or --datasets)", file=sys.stderr)
    plan = [(s, c) for s in steps for c in commands(s, datasets)]
    if a.list or a.dry_run:
        for s, c in plan:
            print(f"{s:14s} oktopy {' '.join(c)}")
        return 0
    if "de-cache" in steps and datasets:
        problem = check_r()
        if problem:
            print(f"oktopy run: {problem}", file=sys.stderr)
            return 1
    t_all = time.time()
    for s, c in plan:
        t0 = time.time()
        print(f"\n=== {s}: oktopy {' '.join(c)}", flush=True)
        r = subprocess.run([sys.executable, "-m", "oktopy", "--project", str(cfg.REPO), *c])
        if r.returncode:
            print(f"oktopy run: step {s} failed (exit {r.returncode}); fix it, then oktopy run --from {s}",
                  file=sys.stderr)
            return r.returncode
        print(f"=== {s}: {time.time() - t0:.0f} s", flush=True)
    print(f"\noktopy run: {len(plan)} commands in {time.time() - t_all:.0f} s; store {cfg.STORE_DB}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
