"""The worked example: a project with two small synthetic datasets and the claims of two fictitious papers, which runs
end to end on a laptop in a few minutes.

Usage: oktopy example DIR [--force] [--no-checks]
then:  cd DIR && oktopy run        (or: snakemake --cores 4)

DIR is made by `oktopy new` and then given the example's field: two regions (front, back), four cell classes with
marker genes, a case-control contrast and a disease stage, two dataset loaders that simulate nuclei with planted
effects (src/example_claims/build/_synthetic.py), and ten claims whose true answer is known. `oktopy run` harmonises
both datasets, builds pseudobulk and the edgeR DE cache, re-tests every claim in both datasets, pools the tests,
rates every claim and writes the store and the atlas (reports/atlas.html). EXAMPLE.md in DIR lists the expected
verdicts.
"""
import argparse
import sys
from pathlib import Path

from oktopy import newproject

HERE = Path(__file__).resolve().parent / "example"
NAME, PACKAGE = "example-claims", "example_claims"
TITLE = "Example claim store (synthetic data)"
AUTHOR = "oktopy worked example"
DESCRIPTION = ("Claims of two fictitious papers on an example disease, re-tested in two synthetic single-nucleus "
               "datasets with planted effects: the oktopy worked example")
# config edits on top of the template (the field's vocabularies, datasets and wording)
EDITS = {
    "config/schema.toml": [('region = ["multiple"]', 'region = ["front", "back", "multiple"]'),
                           ('cell_type = ["All"]', 'cell_type = ["Epithelial", "Immune", "Fibroblast", "Endothelial", "All"]')],
    "config/project.toml": [
        ("datasets = []", 'datasets = ["smith2024", "jones2025"]'),
        ("default_mem_gb = 64", "default_mem_gb = 4"),
        ('examples = { claim = "example2024-C01", finding = "F-example", statement = "CS01", gene = "GENE1", cell = "All", '
         'region = "multiple", dataset = "example2024", background = "BG-example", paper = "example2024", '
         'cell_types = "All", regions = "multiple" }',
         'examples = { claim = "smith2024-C02", finding = "F-example", statement = "CS01", gene = "SPP1", cell = "Immune", '
         'region = "front", dataset = "smith2024", background = "BG-example", paper = "smith2024", '
         'cell_types = "Immune", regions = "front" }'),
        ('keywords = []', 'keywords = ["example", "synthetic data", "claim replication"]'),
        ('disease = ""', 'disease = "example disease"\nexamples = ["Which claims about immune cells hold up?", '
                         '"Is SPP1 up in immune cells?"]'),
    ],
}
REMOVE = ["claims/example2024.yaml", f"src/{PACKAGE}/build/example2024.py"]


def edit(out, rel, pairs):
    p = out / rel
    s = p.read_text()
    for a, b in pairs:
        if a not in s:
            raise SystemExit(f"example: cannot find {a[:60]!r} in {rel} (template changed?)")
        s = s.replace(a, b)
    p.write_text(s)


def make(out, force=False):
    """Write the example project into `out` (a new project plus the example's field); returns its path."""
    values = newproject.values_for(NAME, TITLE, DESCRIPTION, PACKAGE, author=AUTHOR)
    out = newproject.generate(out, values, force)
    for rel in REMOVE:
        (out / rel).unlink()
    newproject.copy_tree(HERE, out, values)
    for rel, pairs in EDITS.items():
        edit(out, rel, pairs)
    newproject.build_schemas(out)  # the example's vocabularies
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(prog="oktopy example", description=__doc__.split("\n")[0])
    ap.add_argument("dir")
    ap.add_argument("--force", action="store_true", help="replace an existing oktopy project in DIR")
    ap.add_argument("--no-checks", action="store_true", help="do not run the project's checks")
    a = ap.parse_args(sys.argv[1:] if argv is None else argv)
    out = make(a.dir, a.force)
    print(f"{out}: the oktopy worked example (two synthetic datasets, two papers, ten claims)")
    ok = None if a.no_checks else newproject.run_checks(out)
    print(f"next: cd {out} && oktopy run      (about 2 minutes; needs oktopy[analysis] and R with edgeR); "
          "EXAMPLE.md lists the expected verdicts")
    return 1 if ok is False else 0


if __name__ == "__main__":
    sys.exit(main())
