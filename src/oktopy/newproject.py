"""Start a new claim-replication project for a field with one command.

Usage: oktopy new DIR --name NAME --title TITLE [--description TEXT] [--package PKG] [--env-prefix PREFIX]
                  [--repository URL] [--author NAME] [--requirement SPEC] [--git] [--force] [--no-checks]

Writes DIR from the project template (placeholders {{name}}, {{title}}, ... filled in), with:
  config/      project.toml (paths, pipeline, release rules, facts about the field), schema.toml (claim format),
               contrasts.toml, harmonise.toml (data schema, cell classes, markers), contexts.toml,
               server_instructions.md
  claims/      an example claim file, empty findings / reference knowledge / questions, SCHEMA.md, RATING.md (the
               engine's rating rules), REGIONS.md, SIGNATURES.md, CURATED_TESTS.md, and the claim schemas built from
               the new schema.toml
  inventory/   papers.tsv and datasets.tsv with an example row
  src/PKG/     a package for the field (--package, default from the name) with an example dataset loader
  workflow/    Snakefile with a laptop profile (default) and a SLURM profile
  tests/       the new project's checks (config, contrasts, schema, claim files), run by its CI workflow
  requirements.txt (oktopy pinned to this version), LICENSE, LICENSE-DATA.md, CITATION.cff, README.md, .gitignore
  .claude/skills/claim-replication/   the method step by step as a Claude Code skill
then runs the new project's checks once (pytest; skipped with a note when pytest is missing). With --git it also makes
the first commit. The engine is not copied: the project uses the installed oktopy, pinned in requirements.txt and
[project] oktopy.
"""
import argparse
import re
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

from oktopy import __version__

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE / "template"
SKILL = HERE / "skill" / "claim-replication"
SKILL_DEST = ".claude/skills/claim-replication"
REPOSITORY = "https://github.com/sergiomarco25/oktopy"
# how a project installs this version (until oktopy is on PyPI: from its repository, at the version's tag)
REQUIREMENT = f"oktopy @ git+{REPOSITORY}@v{__version__}"


def fill(text, values):
    out = re.sub(r"\{\{(\w+)\}\}", lambda m: values.get(m.group(1), m.group(0)), text)
    left = re.findall(r"\{\{\w+\}\}", out)
    if left:
        raise ValueError(f"unfilled placeholders: {sorted(set(left))}")
    return out


def dest_name(rel):
    """Template path -> project path: dot-gitignore -> .gitignore, dot-github/ -> .github/."""
    return "/".join("." + p[4:] if p.startswith("dot-") else p for p in rel.split("/"))


# copied as they are: their double braces are Snakemake / Python syntax, not placeholders
VERBATIM = {"workflow/Snakefile"}


def copy_tree(src, out, values, skip=()):
    """Copy a template tree into `out`, filling placeholders in file names and text files."""
    for f in sorted(p for p in Path(src).rglob("*") if p.is_file() and "__pycache__" not in p.parts):
        rel = f.relative_to(src).as_posix()
        if rel in skip:
            continue
        dst = out / fill(dest_name(rel), values)
        dst.parent.mkdir(parents=True, exist_ok=True)
        if rel in VERBATIM:
            shutil.copy2(f, dst)
            continue
        try:
            dst.write_text(fill(f.read_text(), values))
        except UnicodeDecodeError:
            shutil.copy2(f, dst)


def values_for(name, title, description="", package=None, env_prefix=None, repository="", author="",
               requirement=None):
    if not re.fullmatch(r"[a-z][a-z0-9-]*", name):
        raise SystemExit("--name: lower case letters, digits and '-', starting with a letter")
    package = package or name.replace("-", "_")
    if not re.fullmatch(r"[a-z_][a-z0-9_]*", package) or package in ("oktopy", "claimstore", "tests"):
        raise SystemExit("--package: a Python identifier other than oktopy, claimstore or tests")
    return {"name": name, "title": title.replace('"', "'"), "description": (description or title).replace('"', "'"),
            "package": package, "env_prefix": env_prefix or re.sub(r"\W", "_", name).upper(),
            "repository": repository, "author": author or "", "date": str(date.today()),
            "oktopy_version": __version__, "oktopy_requirement": requirement or REQUIREMENT}


def generate(out, values, force=False):
    out = Path(out).resolve()
    if out.exists() and any(out.iterdir()):
        if not force:
            raise SystemExit(f"{out} exists and is not empty (use --force to replace it)")
        if not (out / "config" / "project.toml").is_file():
            raise SystemExit(f"{out} is not an oktopy project (no config/project.toml); not replacing it")
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    copy_tree(TEMPLATE, out, values)
    shutil.copytree(SKILL, out / SKILL_DEST, ignore=shutil.ignore_patterns("__pycache__"))
    build_schemas(out)
    return out


def build_schemas(out):
    """The claim schemas of the project, built from its own config/schema.toml."""
    r = subprocess.run([sys.executable, "-m", "oktopy", "--project", str(out), "schema"], capture_output=True,
                       text=True)
    if r.returncode:
        raise SystemExit(f"oktopy schema failed in {out}:\n{r.stdout}{r.stderr}")


def run_checks(out):
    missing = [m for m in ("pytest", "yaml", "jsonschema", "pandas") if subprocess.run(
        [sys.executable, "-c", f"import {m}"], capture_output=True).returncode]
    if missing:
        print(f"checks skipped (missing: {', '.join(missing)}); run `pytest tests` in {out}")
        return None
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "tests"], cwd=out,
                       capture_output=True, text=True)
    last = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr[-500:]
    print(f"project checks: {last}")
    if r.returncode:
        print(r.stdout[-3000:])
    return r.returncode == 0


def git_commit(out, message):
    subprocess.run(["git", "init", "-q", "-b", "main"], cwd=out, check=True)
    subprocess.run(["git", "add", "-A"], cwd=out, check=True)
    subprocess.run(["git", "commit", "-q", "-m", message], cwd=out, check=True)
    print("committed")


def parser(prog="oktopy new", doc=__doc__):
    ap = argparse.ArgumentParser(prog=prog, description=doc.split("\n")[0])
    ap.add_argument("dir")
    ap.add_argument("--name", required=True, help="short id, e.g. lung-fibrosis-claims (lower case, digits, '-')")
    ap.add_argument("--title", required=True, help="e.g. 'Lung fibrosis single-cell claim store'")
    ap.add_argument("--description", default="")
    ap.add_argument("--package", help="Python package for the field (default: from the name)")
    ap.add_argument("--env-prefix", help="environment variable prefix (default: from the name, upper case)")
    ap.add_argument("--repository", default="", help="public repository URL, if any")
    ap.add_argument("--author", default="")
    ap.add_argument("--requirement", help=f"how the project installs oktopy (default: {REQUIREMENT})")
    ap.add_argument("--git", action="store_true", help="git init and commit the new project")
    ap.add_argument("--force", action="store_true", help="replace an existing oktopy project in DIR")
    ap.add_argument("--no-checks", action="store_true", help="do not run the new project's checks")
    return ap


def main(argv=None):
    a = parser().parse_args(sys.argv[1:] if argv is None else argv)
    values = values_for(a.name, a.title, a.description, a.package, a.env_prefix, a.repository, a.author,
                        a.requirement)
    out = generate(a.dir, values, a.force)
    print(f"{out}: new project {a.name!r} (package {values['package']}, env prefix {values['env_prefix']}, "
          f"oktopy {__version__})")
    ok = None if a.no_checks else run_checks(out)
    if a.git:
        git_commit(out, f"New claim-replication project (oktopy {__version__})")
    print(f"next: cd {out}; edit config/schema.toml and claims/; see README.md")
    return 1 if ok is False else 0


if __name__ == "__main__":
    sys.exit(main())
