#!/usr/bin/env python3
"""Screen cancer-level multi-route alteration-to-expression support before model development."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

from audit_ucec_expression_signal import hgnc_map
from refine_functional_sources import API, matches


SPECS = {
    "UCEC": {
        "study": "ucec_tcga_pan_can_atlas_2018",
        "genes": {"KRAS": 3845, "PTEN": 5728, "PIK3R1": 5295, "PIK3CA": 5290, "TP53": 7157, "CTNNB1": 1499},
        "routes": [("KRAS_CANONICAL", "MAPK", 1), ("PTEN_LOF", "PI3K", 1), ("PIK3R1_LOF", "PI3K", 1),
                   ("PIK3CA_CANONICAL", "PI3K", 1), ("TP53_DAMAGING", "p53", -1), ("CTNNB1_ACTIVATING", "WNT", 1)],
    },
    "LUAD": {
        "study": "luad_tcga_pan_can_atlas_2018",
        "genes": {"KRAS": 3845, "EGFR": 1956, "TP53": 7157},
        "routes": [("KRAS_CANONICAL", "MAPK", 1), ("EGFR_ACTIVATING", "EGFR", 1),
                   ("EGFR_ACTIVATING", "MAPK", 1), ("TP53_DAMAGING", "p53", -1)],
    },
    "COADREAD": {
        "study": "coadread_tcga_pan_can_atlas_2018",
        "genes": {"APC": 324, "KRAS": 3845, "BRAF": 673, "PIK3CA": 5290, "TP53": 7157},
        "routes": [("APC_LOF", "WNT", 1), ("KRAS_CANONICAL", "MAPK", 1), ("BRAF_V600E", "MAPK", 1),
                   ("PIK3CA_CANONICAL", "PI3K", 1), ("TP53_DAMAGING", "p53", -1)],
    },
    "BRCA": {
        "study": "brca_tcga_pan_can_atlas_2018",
        "genes": {"PIK3CA": 5290, "PTEN": 5728, "AKT1": 207, "TP53": 7157, "ERBB2": 2064},
        "routes": [("PIK3CA_CANONICAL", "PI3K", 1), ("PTEN_LOF", "PI3K", 1), ("AKT1_E17K", "PI3K", 1),
                   ("TP53_DAMAGING", "p53", -1), ("ERBB2_AMP", "EGFR", 1)],
    },
}


def load_progeny(cache: Path) -> pd.DataFrame:
    raw = pd.read_csv(cache, sep="\t")
    wide = raw.pivot_table(index=["genesymbol", "record_id"], columns="label", values="value", aggfunc="first").reset_index()
    wide["p_value"] = pd.to_numeric(wide["p_value"], errors="coerce")
    wide["weight"] = pd.to_numeric(wide["weight"], errors="coerce")
    return wide.dropna(subset=["pathway", "p_value", "weight"])


def extra_match(source: str, gene: str, record: dict) -> bool:
    if source == "CTNNB1_ACTIVATING" and gene == "CTNNB1":
        change = (record.get("proteinChange") or "").replace("p.", "").upper()
        return any(change.startswith(x) for x in ["S33", "S37", "T41", "S45"])
    return matches(source, gene, record)


def clinical(session: requests.Session, study: str, kind: str) -> list[dict]:
    r = session.get(f"{API}/studies/{study}/clinical-data", params={"clinicalDataType": kind, "pageSize": 100000}, timeout=180)
    r.raise_for_status()
    return r.json()


def screen(session: requests.Session, cancer: str, spec: dict, models: dict[str, pd.DataFrame]) -> dict:
    study = spec["study"]
    samples = session.get(f"{API}/studies/{study}/samples", params={"pageSize": 10000}, timeout=90).json()
    primary = [x for x in samples if "Primary" in x.get("sampleType", "")]
    sample_to_patient = {x["sampleId"]: x["patientId"] for x in primary}
    sources = sorted({x[0] for x in spec["routes"]})
    carriers = {x: set() for x in sources}
    mut = session.post(
        f"{API}/molecular-profiles/{study}_mutations/mutations/fetch", params={"projection": "DETAILED"},
        json={"sampleIds": list(sample_to_patient), "entrezGeneIds": list(spec["genes"].values())}, timeout=180,
    )
    mut.raise_for_status()
    for record in mut.json():
        sid = record.get("sampleId")
        if sid not in sample_to_patient:
            continue
        gene = (record.get("gene") or {}).get("hugoGeneSymbol", "")
        for source in sources:
            if extra_match(source, gene, record):
                carriers[source].add(sample_to_patient[sid])
    if "ERBB2_AMP" in carriers:
        cna = session.post(
            f"{API}/molecular-profiles/{study}_gistic/molecular-data/fetch",
            json={"sampleIds": list(sample_to_patient), "entrezGeneIds": [2064]}, timeout=180,
        )
        cna.raise_for_status()
        carriers["ERBB2_AMP"] = {sample_to_patient[x["sampleId"]] for x in cna.json() if x.get("value") is not None and int(float(x["value"])) == 2}

    sc, pc = defaultdict(dict), defaultdict(dict)
    for row in clinical(session, study, "SAMPLE"):
        if row.get("sampleId") in sample_to_patient:
            sc[sample_to_patient[row["sampleId"]]][row["clinicalAttributeId"]] = row.get("value")
    for row in clinical(session, study, "PATIENT"):
        pc[row["patientId"]][row["clinicalAttributeId"]] = row.get("value")

    all_genes = pd.concat(models.values()).drop_duplicates("entrez")
    expr_response = session.post(
        f"{API}/molecular-profiles/{study}_rna_seq_v2_mrna/molecular-data/fetch",
        json={"sampleIds": list(sample_to_patient), "entrezGeneIds": all_genes["entrez"].astype(int).tolist()}, timeout=300,
    )
    expr_response.raise_for_status()
    expr = pd.DataFrame(expr_response.json()).pivot_table(index="sampleId", columns="entrezGeneId", values="value", aggfunc="first")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    z = (expr - expr.mean()) / expr.std(ddof=0)
    scores = {}
    for pathway, model in models.items():
        w = model.drop_duplicates("entrez").set_index("entrez")["weight"].astype(float)
        w = w[w.index.isin(z.columns)]
        scores[pathway] = z[w.index].mul(w, axis=1).sum(axis=1, min_count=max(20, len(w)//2)) / np.sqrt(np.square(w).sum())
    score_frame = pd.DataFrame(scores)
    rows = []
    for sid, values in score_frame.iterrows():
        patient = sample_to_patient[sid]
        s, p = sc[patient], pc[patient]
        rows.append({
            "sample": sid, **values.to_dict(), **{x: int(patient in carriers[x]) for x in sources},
            "SUBTYPE": p.get("SUBTYPE"), "TUMOR_TYPE": s.get("TUMOR_TYPE"),
            "MSI": pd.to_numeric(s.get("MSI_SCORE_MANTIS"), errors="coerce"),
        })
    frame = pd.DataFrame(rows)
    covars = []
    for x in ["SUBTYPE", "TUMOR_TYPE"]:
        if frame[x].notna().mean() >= 0.8 and frame[x].nunique() > 1:
            covars.append(f"C({x})")
    if frame["MSI"].notna().mean() >= 0.8:
        covars.append("MSI")
    results = []
    for source, pathway, expected_sign in spec["routes"]:
        needed = [source, pathway] + [x for x in ["SUBTYPE", "TUMOR_TYPE", "MSI"] if any(x in c for c in covars)]
        use = frame.dropna(subset=needed)
        fit = smf.ols(f"{pathway} ~ {source}" + (" + " + " + ".join(covars) if covars else ""), data=use).fit(cov_type="HC3")
        beta, p2 = float(fit.params[source]), float(fit.pvalues[source])
        pdir = p2 / 2 if beta * expected_sign > 0 else 1 - p2 / 2
        results.append({"source": source, "pathway": pathway, "expected_sign": expected_sign, "carriers": len(carriers[source]),
                        "n": len(use), "beta": beta, "p_two_sided": p2, "p_directional": pdir})
    q = multipletests([x["p_directional"] for x in results], method="fdr_bh")[1]
    for row, value in zip(results, q):
        row["q_directional_within_cancer"] = float(value)
        row["passes"] = bool(row["beta"] * row["expected_sign"] > 0 and value < 0.1 and row["carriers"] >= 20)
    return {"study": study, "samples": len(sample_to_patient), "covariates": covars, "routes": results,
            "supported_routes": sum(x["passes"] for x in results), "supported_pathways": sorted({x["pathway"] for x in results if x["passes"]})}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    mapping = hgnc_map(args.cache_dir / "hgnc_complete_set.txt")
    full = load_progeny(args.cache_dir / "progeny_omnipath.tsv")
    wanted = sorted({route[1] for spec in SPECS.values() for route in spec["routes"]})
    models = {}
    for pathway in wanted:
        part = full[full["pathway"] == pathway].sort_values("p_value").drop_duplicates("genesymbol").head(100).copy()
        part["entrez"] = part["genesymbol"].map(mapping)
        models[pathway] = part.dropna(subset=["entrez"])
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "ARISE-WP0/0.1"})
    results = {c: screen(session, c, s, models) for c, s in SPECS.items()}
    payload = {"results": results, "route_gate": "Expected direction, within-cancer directional FDR<0.10, and >=20 carriers."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
