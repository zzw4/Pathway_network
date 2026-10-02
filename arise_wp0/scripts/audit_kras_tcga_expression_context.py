#!/usr/bin/env python3
"""Test shared and cancer-specific KRAS MAPK transcriptional footprints in TCGA."""

from __future__ import annotations

import argparse
import itertools
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

from audit_ucec_expression_signal import get_progeny_model, hgnc_map
from refine_functional_sources import API, matches


COHORTS = {
    "UCEC": "ucec_tcga_pan_can_atlas_2018",
    "LUAD": "luad_tcga_pan_can_atlas_2018",
    "COADREAD": "coadread_tcga_pan_can_atlas_2018",
}


def clinical(session: requests.Session, study: str, kind: str) -> list[dict]:
    r = session.get(f"{API}/studies/{study}/clinical-data", params={"clinicalDataType": kind, "pageSize": 100000}, timeout=180)
    r.raise_for_status()
    return r.json()


def audit(session: requests.Session, cancer: str, study: str, genes: pd.DataFrame) -> dict:
    samples = session.get(f"{API}/studies/{study}/samples", params={"pageSize": 10000}, timeout=90).json()
    primary = [x for x in samples if "Primary" in x.get("sampleType", "")]
    sample_to_patient = {x["sampleId"]: x["patientId"] for x in primary}
    mut = session.post(
        f"{API}/molecular-profiles/{study}_mutations/mutations/fetch", params={"projection": "DETAILED"},
        json={"sampleIds": list(sample_to_patient), "entrezGeneIds": [3845]}, timeout=180,
    )
    mut.raise_for_status()
    carriers = {x["patientId"] for x in mut.json() if matches("KRAS_CANONICAL", "KRAS", x)}

    sc, pc = defaultdict(dict), defaultdict(dict)
    for row in clinical(session, study, "SAMPLE"):
        if row.get("sampleId") in sample_to_patient:
            sc[sample_to_patient[row["sampleId"]]][row["clinicalAttributeId"]] = row.get("value")
    for row in clinical(session, study, "PATIENT"):
        pc[row["patientId"]][row["clinicalAttributeId"]] = row.get("value")

    profile = f"{study}_rna_seq_v2_mrna"
    r = session.post(
        f"{API}/molecular-profiles/{profile}/molecular-data/fetch",
        json={"sampleIds": list(sample_to_patient), "entrezGeneIds": genes["entrez"].astype(int).tolist()}, timeout=240,
    )
    r.raise_for_status()
    expr = pd.DataFrame(r.json()).pivot_table(index="sampleId", columns="entrezGeneId", values="value", aggfunc="first")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    z = (expr - expr.mean()) / expr.std(ddof=0)
    weights = genes.drop_duplicates("entrez").set_index("entrez")["weight"].astype(float)
    weights = weights[weights.index.isin(z.columns)]
    score = z[weights.index].mul(weights, axis=1).sum(axis=1, min_count=len(weights)//2) / np.sqrt(np.square(weights).sum())

    rows = []
    for sid in z.index:
        patient = sample_to_patient[sid]
        s, p = sc[patient], pc[patient]
        rows.append({
            "sample": sid, "patient": patient, "MAPK": score.get(sid),
            "KRAS_CANONICAL": int(patient in carriers), "SUBTYPE": p.get("SUBTYPE"),
            "TUMOR_TYPE": s.get("TUMOR_TYPE"), "MSI": pd.to_numeric(s.get("MSI_SCORE_MANTIS"), errors="coerce"),
        })
    meta = pd.DataFrame(rows)
    # Keep covariates only when sufficiently complete and nondegenerate.
    terms = ["KRAS_CANONICAL"]
    needed = ["KRAS_CANONICAL"]
    for cov in ["SUBTYPE", "TUMOR_TYPE"]:
        if meta[cov].notna().mean() >= 0.8 and meta[cov].nunique() > 1:
            terms.append(f"C({cov})")
            needed.append(cov)
    if meta["MSI"].notna().mean() >= 0.8:
        terms.append("MSI")
        needed.append("MSI")
    use = meta.dropna(subset=needed + ["MAPK"]).copy()
    module_fit = smf.ols(f"MAPK ~ {' + '.join(terms)}", data=use).fit(cov_type="HC3")

    gene_results = []
    symbol_by_entrez = genes.drop_duplicates("entrez").set_index("entrez")["genesymbol"].to_dict()
    for entrez in weights.index:
        outcome = z[entrez].rename("Y")
        frame = use.merge(outcome, left_on="sample", right_index=True).dropna(subset=["Y"])
        fit = smf.ols(f"Y ~ {' + '.join(terms)}", data=frame).fit(cov_type="HC3")
        gene_results.append({
            "gene": symbol_by_entrez[entrez], "entrez": int(entrez), "n": int(len(frame)),
            "beta": float(fit.params["KRAS_CANONICAL"]), "p": float(fit.pvalues["KRAS_CANONICAL"]),
        })
    qvals = multipletests([x["p"] for x in gene_results], method="fdr_bh")[1]
    for row, q in zip(gene_results, qvals):
        row["q"] = float(q)
    return {
        "study": study, "primary_samples": len(sample_to_patient), "complete_cases": len(use),
        "kras_carriers": len(carriers), "mapk_genes": len(weights), "covariates": terms,
        "module": {
            "beta": float(module_fit.params["KRAS_CANONICAL"]),
            "ci95": [float(x) for x in module_fit.conf_int().loc["KRAS_CANONICAL"]],
            "p": float(module_fit.pvalues["KRAS_CANONICAL"]),
        },
        "gene_results": gene_results,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    model = get_progeny_model(args.cache_dir / "progeny_omnipath.tsv")
    mapping = hgnc_map(args.cache_dir / "hgnc_complete_set.txt")
    genes = model[model["pathway"] == "MAPK"].sort_values("p_value").drop_duplicates("genesymbol").head(100).copy()
    genes["entrez"] = genes["genesymbol"].map(mapping)
    genes = genes.dropna(subset=["entrez"])
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "ARISE-WP0/0.1"})
    results = {c: audit(session, c, s, genes) for c, s in COHORTS.items()}

    comparisons = {}
    for a, b in itertools.combinations(results, 2):
        left = pd.DataFrame(results[a]["gene_results"])[["gene", "beta", "q"]].rename(columns={"beta": "a_beta", "q": "a_q"})
        right = pd.DataFrame(results[b]["gene_results"])[["gene", "beta", "q"]].rename(columns={"beta": "b_beta", "q": "b_q"})
        merged = left.merge(right, on="gene")
        comparisons[f"{a}_vs_{b}"] = {
            "genes": len(merged), "beta_pearson": float(merged["a_beta"].corr(merged["b_beta"])),
            "same_direction": int((np.sign(merged["a_beta"]) == np.sign(merged["b_beta"])).sum()),
            "significant_both_q05": int(((merged["a_q"] < 0.05) & (merged["b_q"] < 0.05)).sum()),
            "significant_either_q05": int(((merged["a_q"] < 0.05) | (merged["b_q"] < 0.05)).sum()),
        }
    payload = {"results": results, "pairwise_context_comparisons": comparisons}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    compact = {c: {k: v for k, v in x.items() if k != "gene_results"} for c, x in results.items()}
    print(json.dumps({"cohorts": compact, "comparisons": comparisons}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
