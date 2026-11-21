"""Synthetic single-nucleus counts for the oktopy worked example (no download; the same numbers on every machine).

Each dataset has donors (case or control), one or two regions, four cell classes with marker genes, and effects
planted in the case donors, so every claim of the example has a known answer:

  Immune share             x2 in cases (epithelial x0.5)     claims smith2024-C01, jones2025-C01  (true)
  SPP1 in Immune           x3 in cases                       smith2024-C02                         (true)
  interferon genes         x2 (smith2024), x1.6 (jones2025)  smith2024-C03, jones2025-C03          (true)
    (IFI6, ISG15, MX1, OAS1, IFIT3) in Epithelial
  HIF1A in Endothelial     x0.5 in cases                     smith2024-C04 claims up               (false: down)
  Endothelial share        unchanged (centred log-ratio)     smith2024-C05 claims no difference    (true)
  MKI67 in Epithelial      x1.4 per stage step (stage 0-3)   smith2024-C06 claims a rise with stage (true)
  COL3A1 in Fibroblast     unchanged                         jones2025-C02 claims down             (false: no change)

Counts are negative binomial (gamma-Poisson) around a class profile, with donor and region noise.
"""
import numpy as np
import pandas as pd

MARKERS = {"Epithelial": ["EPCAM", "KRT8", "KRT18", "CDH1", "KRT19"],
           "Immune": ["PTPRC", "CD74", "LYZ", "CD3E", "CD68"],
           "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"],
           "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "FLT1"]}
AUTHOR_LABEL = {"Epithelial": "Epi", "Immune": "Imm", "Fibroblast": "Fib", "Endothelial": "Endo"}
AUTHOR_TO_CLASS = {v: k for k, v in AUTHOR_LABEL.items()}
IFN = ["IFI6", "ISG15", "MX1", "OAS1", "IFIT3"]
EFFECT_GENES = ["SPP1", *IFN, "HIF1A", "COL3A1", "MKI67"]
N_FILLER = 600
SHARES = {"Epithelial": 0.45, "Immune": 0.15, "Fibroblast": 0.28, "Endothelial": 0.12}


def genes():
    markers = [g for gs in MARKERS.values() for g in gs]
    return markers + EFFECT_GENES + [f"GENE{i:04d}" for i in range(1, N_FILLER + 1)]


def make(dataset, regions, n_case, n_control, seed, ifn_fold=2.0, cells_per_unit=150, umis=2500):
    """AnnData with raw counts and the obs columns of the harmonised schema."""
    import anndata as ad
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    g = genes()
    gi = {x: i for i, x in enumerate(g)}
    classes = list(MARKERS)
    # class expression profiles (log scale): shared baseline, class-specific noise, markers high in their class
    base = rng.normal(-1.0, 1.1, len(g))
    for x in EFFECT_GENES:
        base[gi[x]] = rng.normal(0.8, 0.3)
    prof = {c: base + rng.normal(0, 0.35, len(g)) for c in classes}
    for c, ms in MARKERS.items():
        for other in classes:
            for m in ms:
                prof[other][gi[m]] = base[gi[m]] + (3.0 if other == c else -1.5)

    donors = [f"{dataset}_D{i + 1:02d}" for i in range(n_case + n_control)]
    disease = ["case"] * n_case + ["control"] * n_control
    meta = pd.DataFrame({
        "donor": donors, "disease": disease,
        "disease_raw": ["Example disease" if d == "case" else "Healthy" for d in disease],
        "sex": rng.choice(["F", "M"], len(donors)),
        "age": rng.integers(45, 86, len(donors)).astype(float),
        "pmi": np.round(rng.uniform(2, 24, len(donors)), 1),
        "stage": [float(rng.integers(1, 4)) if d == "case" else 0.0 for d in disease]})
    region_eff = {r: rng.normal(0, 0.15, len(g)) for r in regions}

    blocks, obs = [], []
    for d in meta.itertuples():
        donor_eff = rng.normal(0, 0.15, len(g))
        # cases: immune cells x2 at the expense of epithelial cells (x0.5), so the geometric mean of the shares and
        # with it the centred log-ratio of fibroblasts and endothelial cells stay as in controls
        fold = {"Immune": 2.0, "Epithelial": 0.5} if d.disease == "case" else {}
        shares = np.array([SHARES[c] * fold.get(c, 1.0) for c in classes])
        for reg in regions:
            p = rng.dirichlet(shares / shares.sum() * 80)
            n = rng.multinomial(rng.poisson(cells_per_unit), p)
            for c, k in zip(classes, n):
                if not k:
                    continue
                logmu = prof[c] + donor_eff + region_eff[reg]
                if d.disease == "case":
                    if c == "Immune":
                        logmu[gi["SPP1"]] += np.log(3.0)
                    if c == "Epithelial":
                        logmu[[gi[x] for x in IFN]] += np.log(ifn_fold)
                    if c == "Endothelial":
                        logmu[gi["HIF1A"]] -= np.log(2.0)
                if c == "Epithelial":
                    logmu[gi["MKI67"]] += np.log(1.4) * d.stage
                mu = np.exp(logmu)
                mu = mu / mu.sum()
                lib = rng.lognormal(np.log(umis), 0.3, k)[:, None] * mu[None, :]
                lam = rng.gamma(1 / 0.3, lib * 0.3)  # negative binomial, dispersion 0.3
                blocks.append(sp.csr_matrix(rng.poisson(lam).astype(np.float32)))
                obs.append(pd.DataFrame({"donor": d.donor, "sample": f"{d.donor}_{reg}", "region": reg,
                                         "author_cell_type": AUTHOR_LABEL[c]}, index=range(k)))
    X = sp.vstack(blocks).tocsr()
    o = pd.concat(obs, ignore_index=True).merge(meta, on="donor", how="left")
    o.index = [f"{dataset}_c{i}" for i in range(len(o))]
    o["dataset"] = dataset
    o["cell_class"] = o.author_cell_type.map(AUTHOR_TO_CLASS)
    return ad.AnnData(X, obs=o, var=pd.DataFrame(index=g))
