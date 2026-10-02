#!/usr/bin/env python3
"""ARISE Gate B: source-specific observed-vs-prior direction audit in TCGA-UCEC."""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import statsmodels.api as sm
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

from audit_ucec_expression_signal import clinical, hgnc_map
from refine_functional_sources import API, matches


STUDY = "ucec_tcga_pan_can_atlas_2018"
PROFILE = f"{STUDY}_rna_seq_v2_mrna"
GENES = {"KRAS": 3845, "CTNNB1": 1499, "PTEN": 5728, "PIK3CA": 5290, "PIK3R1": 5295, "TP53": 7157}
CTNNB1_RESIDUES = {32, 33, 34, 37, 41, 45}


def ctnnb1_activating(record: dict) -> bool:
    change = (record.get("proteinChange") or "").replace("p.", "").upper()
    match = re.match(r"^[A-Z](\d+)", change)
    return bool(match and int(match.group(1)) in CTNNB1_RESIDUES)


def design_matrix(meta: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    base = meta[["KRAS", "CTNNB1", "PTEN_LOF", "PIK3CA_CANONICAL", "PIK3R1_LOF", "TP53_DAMAGING"]].astype(float)
    used = list(base.columns)
    if meta["MSI"].notna().mean() >= 0.8 and meta["MSI"].nunique(dropna=True) > 1:
        base["MSI"] = pd.to_numeric(meta["MSI"], errors="coerce")
        used.append("MSI")
    for covariate in ("SUBTYPE", "TUMOR_TYPE", "GRADE"):
        if meta[covariate].notna().mean() >= 0.8 and meta[covariate].nunique(dropna=True) > 1:
            dummies = pd.get_dummies(meta[covariate], prefix=covariate, drop_first=True, dtype=float)
            base = pd.concat([base, dummies], axis=1)
            used.append(f"categorical:{covariate}")
    base = sm.add_constant(base, has_constant="add")
    return base, used


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--cache-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    mapping = hgnc_map(args.cache_dir / "hgnc_complete_set.txt")
    expected = {}
    for source in ("KRAS", "CTNNB1"):
        expected[source] = {
            row["target"]: int(row["predicted_sign"])
            for row in manifest["sources"][source]["directional_targets"]
            if row["target"] != source and row["target"] in mapping
        }
    requested_symbols = sorted(set().union(*[set(x) for x in expected.values()]))
    requested_entrez = sorted({mapping[x] for x in requested_symbols})
    symbol_by_entrez = {mapping[x]: x for x in requested_symbols}

    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "ARISE-WP0/0.2"})
    samples_response = session.get(f"{API}/studies/{STUDY}/samples", params={"pageSize": 10000}, timeout=90)
    samples_response.raise_for_status()
    primary = [x for x in samples_response.json() if "Primary" in x.get("sampleType", "")]
    sample_to_patient = {x["sampleId"]: x["patientId"] for x in primary}

    mut_response = session.post(
        f"{API}/molecular-profiles/{STUDY}_mutations/mutations/fetch",
        params={"projection": "DETAILED"},
        json={"sampleIds": list(sample_to_patient), "entrezGeneIds": list(GENES.values())}, timeout=240,
    )
    mut_response.raise_for_status()
    carriers = {name: set() for name in ("KRAS", "CTNNB1", "PTEN_LOF", "PIK3CA_CANONICAL", "PIK3R1_LOF", "TP53_DAMAGING")}
    for record in mut_response.json():
        sid = record.get("sampleId")
        if sid not in sample_to_patient:
            continue
        patient = sample_to_patient[sid]
        gene = (record.get("gene") or {}).get("hugoGeneSymbol", "")
        if gene == "KRAS" and matches("KRAS_CANONICAL", gene, record):
            carriers["KRAS"].add(patient)
        if gene == "CTNNB1" and ctnnb1_activating(record):
            carriers["CTNNB1"].add(patient)
        for label in ("PTEN_LOF", "PIK3CA_CANONICAL", "PIK3R1_LOF", "TP53_DAMAGING"):
            if matches(label, gene, record):
                carriers[label].add(patient)

    sample_clinical, patient_clinical = defaultdict(dict), defaultdict(dict)
    for row in clinical(session, "SAMPLE"):
        if row.get("sampleId") in sample_to_patient:
            sample_clinical[sample_to_patient[row["sampleId"]]][row["clinicalAttributeId"]] = row.get("value")
    for row in clinical(session, "PATIENT"):
        patient_clinical[row["patientId"]][row["clinicalAttributeId"]] = row.get("value")

    expr_response = session.post(
        f"{API}/molecular-profiles/{PROFILE}/molecular-data/fetch",
        json={"sampleIds": list(sample_to_patient), "entrezGeneIds": requested_entrez}, timeout=300,
    )
    expr_response.raise_for_status()
    expr_long = pd.DataFrame(expr_response.json())
    expression = expr_long.pivot_table(index="sampleId", columns="entrezGeneId", values="value", aggfunc="first")
    expression = expression.apply(pd.to_numeric, errors="coerce")
    expression = (expression - expression.mean()) / expression.std(ddof=0)

    rows = []
    for sid in expression.index:
        patient = sample_to_patient[sid]
        sc, pc = sample_clinical[patient], patient_clinical[patient]
        rows.append({
            "sample": sid, "patient": patient,
            **{label: int(patient in members) for label, members in carriers.items()},
            "SUBTYPE": pc.get("SUBTYPE"), "TUMOR_TYPE": sc.get("TUMOR_TYPE"), "GRADE": sc.get("GRADE"),
            "MSI": pd.to_numeric(sc.get("MSI_SCORE_MANTIS"), errors="coerce"),
        })
    meta = pd.DataFrame(rows).set_index("sample")
    design, covariates = design_matrix(meta)

    results = {}
    all_gene_rows = []
    for source in ("KRAS", "CTNNB1"):
        gene_rows = []
        for symbol, expected_sign in expected[source].items():
            entrez = mapping[symbol]
            if entrez not in expression.columns:
                continue
            joined = design.join(expression[entrez].rename("Y"), how="inner").dropna()
            if joined[source].sum() < 10 or joined[source].nunique() < 2:
                continue
            fit = sm.OLS(joined["Y"], joined.drop(columns="Y")).fit(cov_type="HC3")
            beta, p = float(fit.params[source]), float(fit.pvalues[source])
            gene_rows.append({
                "source": source, "gene": symbol, "entrez": int(entrez), "n": int(len(joined)),
                "carriers": int(joined[source].sum()), "expected_sign": expected_sign,
                "beta": beta, "p": p, "direction_match": bool(np.sign(beta) == expected_sign),
            })
        qvals = multipletests([x["p"] for x in gene_rows], method="fdr_bh")[1] if gene_rows else []
        for row, q in zip(gene_rows, qvals):
            row["q"] = float(q)
        all_gene_rows.extend(gene_rows)

        def summary(subset):
            count = len(subset)
            matched = sum(x["direction_match"] for x in subset)
            return {
                "genes": count, "matched": int(matched),
                "concordance": float(matched / count) if count else None,
                "binomial_p_vs_0.5": float(binomtest(matched, count, 0.5, alternative="greater").pvalue) if count else None,
            }

        results[source] = {
            "manifest_targets_mapped": len(expected[source]),
            "tested_targets": len(gene_rows),
            "carriers_in_expression_set": int(meta[source].sum()),
            "all": summary(gene_rows),
            "nominal_p_lt_0.05": summary([x for x in gene_rows if x["p"] < 0.05]),
            "fdr_q_lt_0.10": summary([x for x in gene_rows if x["q"] < 0.10]),
            "median_signed_effect": float(np.median([x["beta"] * x["expected_sign"] for x in gene_rows])) if gene_rows else None,
        }

    payload = {
        "purpose": "ARISE Gate B observed direction audit; sources are never pooled",
        "study": STUDY,
        "model": "per-gene OLS with HC3 robust covariance",
        "covariates": covariates,
        "carrier_counts": {label: int(meta[label].sum()) for label in carriers},
        "results": results,
        "interpretation_policy": {
            "primary_readout": "effect direction concordance on preregistered source-reachable direct TF targets",
            "not_claimed": "causality or alteration-removal attribution",
            "next_gate": "Gate C additive source-specific baselines only if Gate B is non-null and biologically coherent",
        },
        "gene_results": all_gene_rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"carrier_counts": payload["carrier_counts"], "results": results}, indent=2))


if __name__ == "__main__":
    main()
