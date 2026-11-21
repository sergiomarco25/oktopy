"""Method registry: which methods the engine runs, or could run, for each analysis task, ranked into tiers
(oktopy/methods.toml; ranking criteria and rules in its header).

  oktopy registry            validate and print a summary (exit 1 on errors)
  oktopy registry --md OUT   also write the registry as a Markdown page
  oktopy registry --json     the registry with the project's switches applied, as JSON

A project switches methods off in config/methods.toml:
  [disabled]
  ids = ["nebula"]          # never a primary; a project cannot promote a method
Standard library only.
"""
import json
import re
import sys
import tomllib
from pathlib import Path

from oktopy import config as cfg

HERE = Path(__file__).resolve().parent
REGISTRY_FILE = HERE / "methods.toml"
PROJECT_FILE = cfg.REPO / "config" / "methods.toml"
SCHEMA_FILE = cfg.REPO / "config" / "schema.toml"
TIERS = ("primary", "sensitivity", "exploratory", "excluded")
STATUS = ("implemented", "planned")
TASK_KEYS = {"title", "modality", "book", "confidence", "note", "method"}
METHOD_KEYS = {"id", "tier", "status", "engine", "serves", "unit", "requires", "evidence", "why"}
REQUIRED = ("id", "tier", "status", "serves", "unit", "evidence", "why")


def _read(path):
    with open(path, "rb") as fh:
        return tomllib.load(fh)


def load(path=REGISTRY_FILE, project_file=PROJECT_FILE):
    """The registry with the project's disabled methods marked (`disabled = True`), never removed."""
    reg = _read(path)
    off = set()
    if project_file and Path(project_file).is_file():
        off = set(_read(project_file).get("disabled", {}).get("ids", []))
    for task in reg.get("task", {}).values():
        for m in task.get("method", []):
            m["disabled"] = m["id"] in off
    reg["disabled_ids"] = sorted(off)
    return reg


def templates():
    """Claim test templates of the project (config/schema.toml [vocab] test_template), or None without a schema."""
    if not SCHEMA_FILE.is_file():
        return None
    return set(_read(SCHEMA_FILE).get("vocab", {}).get("test_template", []))


def validate(reg, known_templates=None):
    """List of problems (empty = valid)."""
    errs = []
    if not re.fullmatch(r"\d+\.\d+\.\d+", str(reg.get("version", ""))):
        errs.append("version: x.y.z required")
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", str(reg.get("reviewed", ""))):
        errs.append("reviewed: YYYY-MM-DD required")
    if not reg.get("book"):
        errs.append("book: the best-practices release the registry follows")
    tasks = reg.get("task", {})
    if not tasks:
        errs.append("no tasks")
    for name, t in tasks.items():
        extra = set(t) - TASK_KEYS
        if extra:
            errs.append(f"{name}: unknown keys {sorted(extra)}")
        for k in ("title", "modality", "confidence"):
            if not t.get(k):
                errs.append(f"{name}: {k} missing")
        ms = t.get("method", [])
        if not ms:
            errs.append(f"{name}: no methods")
        ids = [m.get("id") for m in ms]
        if len(ids) != len(set(ids)):
            errs.append(f"{name}: duplicate method ids")
        prim = [m for m in ms if m.get("tier") == "primary"]
        if len(prim) > 1:
            errs.append(f"{name}: {len(prim)} primaries (one at most)")
        for m in prim:
            if m.get("disabled"):
                errs.append(f"{name}/{m['id']}: a project cannot switch off the primary")
        for m in ms:
            where = f"{name}/{m.get('id')}"
            extra = set(m) - METHOD_KEYS - {"disabled"}
            if extra:
                errs.append(f"{where}: unknown keys {sorted(extra)}")
            for k in REQUIRED:
                if k not in m or (k != "serves" and not m[k]):
                    errs.append(f"{where}: {k} missing")
            if m.get("tier") not in TIERS:
                errs.append(f"{where}: tier must be one of {TIERS}")
            if m.get("status") not in STATUS:
                errs.append(f"{where}: status must be one of {STATUS}")
            if m.get("status") == "implemented" and not m.get("engine"):
                errs.append(f"{where}: implemented methods name their engine code")
            if known_templates is not None:
                bad = [s for s in m.get("serves", []) if s not in known_templates]
                if bad:
                    errs.append(f"{where}: serves unknown templates {bad}")
    if known_templates is not None:
        served = {s for t in tasks.values() for m in t.get("method", [])
                  if m.get("tier") == "primary" and m.get("status") == "implemented" for s in m.get("serves", [])}
        # every template the engine tests today must have an implemented primary somewhere
        engine_templates = {"composition_donor", "de_pseudobulk", "geneset_shift", "state_signature",
                            "genetic_enrichment", "pathology_correlation", "trajectory", "communication_liana"}
        missing = sorted((engine_templates & known_templates) - served)
        if missing:
            errs.append(f"templates tested today without an implemented primary: {missing}")
    return errs


def rows(reg):
    """One dict per task x method, for tables."""
    out = []
    for name, t in reg.get("task", {}).items():
        for m in t.get("method", []):
            out.append({"task": name, "task_title": t.get("title", ""), "modality": t.get("modality", ""),
                        "confidence": t.get("confidence", ""), "book": t.get("book", ""), **m})
    return out


def markdown(reg):
    order = {k: i for i, k in enumerate(TIERS)}
    lines = [f"# Method registry {reg['version']}", "",
             f"Reviewed {reg['reviewed']}. Follows {reg['book']} ({reg.get('book_url', '')}) and independent benchmarks. "
             "Primary methods set the verdict; sensitivity methods can only lower certainty; exploratory methods are "
             "reported, not rated; excluded methods are not run (except as the original-method diagnostic).", ""]
    if reg.get("disabled_ids"):
        lines += [f"Switched off in this project: {', '.join(reg['disabled_ids'])}.", ""]
    for name, t in reg.get("task", {}).items():
        lines += [f"## {t['title']} (`{name}`)", "",
                  f"Modality: {t['modality']}. Book chapter: {t.get('book') or 'not covered'}. "
                  f"Confidence in the ranking: {t['confidence']}.", ""]
        if t.get("note"):
            lines += [t["note"], ""]
        lines += ["| tier | method | status | unit | why | evidence |", "|---|---|---|---|---|---|"]
        for m in sorted(t.get("method", []), key=lambda m: order.get(m["tier"], 9)):
            st = m["status"] + (", switched off here" if m.get("disabled") else "")
            lines.append(f"| {m['tier']} | {m['id']} | {st} | {m['unit']} | {m['why']} | {'; '.join(m['evidence'])} |")
        lines.append("")
    if reg.get("watch"):
        lines += ["## Sources watched for changes", ""] + [f"- {w}" for w in reg["watch"]] + [""]
    return "\n".join(lines)


def summary(reg):
    rs = rows(reg)
    by = {}
    for r in rs:
        by.setdefault((r["tier"], r["status"]), 0)
        by[(r["tier"], r["status"])] += 1
    parts = [f"{t}: " + ", ".join(f"{by.get((t, s), 0)} {s}" for s in STATUS) for t in TIERS]
    return (f"registry {reg['version']} ({reg['book']}): {len(reg.get('task', {}))} tasks, {len(rs)} methods\n  "
            + "\n  ".join(parts))


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    reg = load()
    errs = validate(reg, templates())
    if "--json" in argv:
        print(json.dumps(reg, indent=1))
    else:
        print(summary(reg))
    if "--md" in argv:
        out = Path(argv[argv.index("--md") + 1])
        out.write_text(markdown(reg))
        print("wrote", out)
    for e in errs:
        print("ERROR", e, file=sys.stderr)
    return 1 if errs else 0


if __name__ == "__main__":
    sys.exit(main())
