#!/usr/bin/env python3
"""Preregistered directional CPTAC-UCEC phosphosite gate for PIK3R1 and KRAS sources."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import statsmodels.formula.api as smf
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

from refine_functional_sources import API, matches


STUDY = "ucec_cptac_2020"
MUTATION_PROFILE = f"{STUDY}_mutations"
PROTEIN_PROFILE = f"{STUDY}_protein_quantification"
PHOSPHO_PROFILE = f"{STUDY}_phosphoproteome"
GENES = {
    "PIK3R1": 5295, "PIK3CA": 5290, "PTEN": 5728, "KRAS": 3845, "AKT2": 208, "GSK3B": 2932,
    "FOXO3": 2309, "TSC2": 7249, "EIF4EBP1": 1978, "RPS6": 6194,
    "MAPK1": 5594, "MAPK3": 5595,
}

# Frozen before association testing. All signs mean higher phosphorylation is expected.
SITE_PANEL = {
    "PI3K": [
        "AKT2_S474", "GSK3B_S9", "FOXO3_T32", "TSC2_T1462",
        "EIF4EBP1_T37", "EIF4EBP1_T46", "EIF4EBP1_S65", "EIF4EBP1_T70",
        "RPS6_S235", "RPS6_S236", "RPS6_S240", "RPS6_S244",
    ],
    "MAPK": ["MAPK1_T185", "MAPK1_Y187", "MAPK3_T202", "MAPK3_Y204"],
}
SOURCE_FOR_MODULE = {"PI3K": "PIK3R1_LOF", "MAPK": "KRAS_CANONICAL"}


def fetch_clinical(session: requests.Session, kind: str) -> list[dict]:
    r = session.get(
        f"{API}/studies/{STUDY}/clinical-data",
        params={"clinicalDataType": kind, "pageSize": 100000}, timeout=180,
    )
    r.raise_for_status()
    return r.json()


def site_gene(site: str) -> str:
    return site.rsplit("_", 1)[0]


def robust_fit(data: pd.DataFrame, outcome: str, source: str, alternate: str, include_protein: bool) -> dict:
    terms = [source, alternate]
    if include_protein:
        terms.append("TOTAL_PROTEIN")
    if data["GENOMICS_SUBTYPE"].notna().sum() >= 0.8 * len(data):
        terms.append("C(GENOMICS_SUBTYPE)")
    if data["PURITY_CANCER"].notna().sum() >= 0.8 * len(data):
        terms.append("PURITY_CANCER")
    needed = [outcome, source, alternate] + (["TOTAL_PROTEIN"] if include_protein else [])
    needed += [x for x in ["GENOMICS_SUBTYPE", "PURITY_CANCER"] if any(x in t for t in terms)]
    use = data.dropna(subset=needed).copy()
    if use[source].sum() < 5 or (len(use) - use[source].sum()) < 5:
        return {"n": int(len(use)), "carriers": int(use[source].sum()), "error": "insufficient source groups"}
    model = smf.ols(f"{outcome} ~ {' + '.join(terms)}", data=use).fit(cov_type="HC3")
    beta = float(model.params[source])
    p_two = float(model.pvalues[source])
    p_directional = p_two / 2 if beta > 0 else 1 - p_two / 2
    return {
        "n": int(len(use)), "carriers": int(use[source].sum()), "beta": beta,
        "se": float(model.bse[source]), "p_two_sided": p_two,
        "p_directional_positive": float(p_directional),
        "ci95": [float(x) for x in model.conf_int().loc[source]],
        "covariates": terms,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "ARISE-WP0/0.1"})

    samples = session.get(f"{API}/studies/{STUDY}/samples", params={"pageSize": 10000}, timeout=90).json()
    sample_ids = [x["sampleId"] for x in samples]
    sample_to_patient = {x["sampleId"]: x["patientId"] for x in samples}

    mutation_response = session.post(
        f"{API}/molecular-profiles/{MUTATION_PROFILE}/mutations/fetch",
        params={"projection": "DETAILED"},
        json={"sampleIds": sample_ids, "entrezGeneIds": [GENES[x] for x in ["PIK3R1", "PIK3CA", "PTEN", "KRAS"]]}, timeout=180,
    )
    mutation_response.raise_for_status()
    carriers = {"PIK3R1_LOF": set(), "PIK3CA_CANONICAL": set(), "PTEN_LOF": set(), "KRAS_CANONICAL": set()}
    for record in mutation_response.json():
        gene = (record.get("gene") or {}).get("hugoGeneSymbol", "")
        for source in carriers:
            if matches(source, gene, record):
                carriers[source].add(record["patientId"])

    sample_clin, patient_clin = defaultdict(dict), defaultdict(dict)
    for row in fetch_clinical(session, "SAMPLE"):
        sample_clin[row["sampleId"]][row["clinicalAttributeId"]] = row.get("value")
    for row in fetch_clinical(session, "PATIENT"):
        patient_clin[row["patientId"]][row["clinicalAttributeId"]] = row.get("value")

    phospho = {}
    for sites in SITE_PANEL.values():
        for site in sites:
            r = session.get(f"{API}/generic-assay-data/{PHOSPHO_PROFILE}/generic-assay/{site}", timeout=120)
            r.raise_for_status()
            phospho[site] = {x["sampleId"]: pd.to_numeric(x.get("value"), errors="coerce") for x in r.json()}

    protein_genes = sorted({site_gene(s) for sites in SITE_PANEL.values() for s in sites})
    protein_response = session.post(
        f"{API}/molecular-profiles/{PROTEIN_PROFILE}/molecular-data/fetch",
        params={"projection": "DETAILED"},
        json={"sampleIds": sample_ids, "entrezGeneIds": [GENES[g] for g in protein_genes]}, timeout=180,
    )
    protein_response.raise_for_status()
    protein = defaultdict(dict)
    for row in protein_response.json():
        gene = (row.get("gene") or {}).get("hugoGeneSymbol", "")
        protein[gene][row["sampleId"]] = pd.to_numeric(row.get("value"), errors="coerce")

    base_rows = []
    for sid in sample_ids:
        patient = sample_to_patient[sid]
        sc, pc = sample_clin[sid], patient_clin[patient]
        base_rows.append({
            "sample": sid, "patient": patient,
            **{source: int(patient in members) for source, members in carriers.items()},
            "GENOMICS_SUBTYPE": sc.get("GENOMICS_SUBTYPE") or pc.get("GENOMICS_SUBTYPE"),
            "PURITY_CANCER": pd.to_numeric(sc.get("PURITY_CANCER") or pc.get("PURITY_CANCER"), errors="coerce"),
        })
    base = pd.DataFrame(base_rows)

    site_results = []
    residual_columns = {"PI3K": [], "MAPK": []}
    module_frame = base.copy()
    for module, sites in SITE_PANEL.items():
        source = SOURCE_FOR_MODULE[module]
        alternate = "KRAS_CANONICAL" if source == "PIK3R1_LOF" else "PIK3R1_LOF"
        for idx, site in enumerate(sites):
            gene = site_gene(site)
            frame = base.copy()
            frame["PHOSPHO"] = frame["sample"].map(phospho[site])
            frame["TOTAL_PROTEIN"] = frame["sample"].map(protein[gene])
            fit = robust_fit(frame, "PHOSPHO", source, alternate, include_protein=True)
            fit.update({"module": module, "site": site, "gene": gene, "source": source})
            site_results.append(fit)

            # Protein-adjusted standardized site for the preregistered module summary.
            available = frame.dropna(subset=["PHOSPHO", "TOTAL_PROTEIN"])
            if len(available) >= 30:
                residual_model = smf.ols("PHOSPHO ~ TOTAL_PROTEIN", data=available).fit()
                resid = pd.Series(residual_model.resid.to_numpy(), index=available["sample"].to_numpy())
                z = (resid - resid.mean()) / resid.std(ddof=0)
                col = f"{module}_{idx}"
                module_frame[col] = module_frame["sample"].map(z)
                residual_columns[module].append(col)

    valid_indices = [i for i, x in enumerate(site_results) if "p_directional_positive" in x]
    if valid_indices:
        qvals = multipletests([site_results[i]["p_directional_positive"] for i in valid_indices], method="fdr_bh")[1]
        for i, q in zip(valid_indices, qvals):
            site_results[i]["q_directional_panel"] = float(q)

    module_results = {}
    for module, columns in residual_columns.items():
        module_frame[f"{module}_MODULE"] = module_frame[columns].mean(axis=1, skipna=True)
        module_frame.loc[module_frame[columns].notna().sum(axis=1) < max(2, len(columns) // 2), f"{module}_MODULE"] = np.nan
        source = SOURCE_FOR_MODULE[module]
        alternate = "KRAS_CANONICAL" if source == "PIK3R1_LOF" else "PIK3R1_LOF"
        module_results[module] = robust_fit(module_frame, f"{module}_MODULE", source, alternate, include_protein=False)
        module_results[module]["source"] = source
        module_results[module]["sites_contributing"] = columns

    exploratory_pi3k_source_screen = {}
    for source in ["PIK3R1_LOF", "PIK3CA_CANONICAL", "PTEN_LOF"]:
        exploratory_pi3k_source_screen[source] = robust_fit(
            module_frame, "PI3K_MODULE", source, "KRAS_CANONICAL", include_protein=False
        )

    expected_positive = [x for x in site_results if x.get("beta", 0) > 0]
    payload = {
        "study": STUDY, "samples": len(sample_ids),
        "carrier_counts": {k: len(v) for k, v in carriers.items()},
        "clinical_completeness": {
            "GENOMICS_SUBTYPE": int(base["GENOMICS_SUBTYPE"].notna().sum()),
            "PURITY_CANCER": int(base["PURITY_CANCER"].notna().sum()),
        },
        "frozen_site_panel": SITE_PANEL, "site_results": site_results, "module_results": module_results,
        "exploratory_pi3k_source_screen": exploratory_pi3k_source_screen,
        "panel_summary": {
            "directionally_positive_sites": len(expected_positive), "tested_sites": len(valid_indices),
            "directional_fdr_lt_0.1": sum(x.get("q_directional_panel", 1) < 0.1 for x in site_results),
        },
        "gate_rule": "Advance only if each source has a positive adjusted module effect and at least one biologically matching site with directional FDR < 0.10; otherwise revise the source or cancer before neural-network development.",
        "limitations": [
            "Observational phosphoproteomics does not establish causality.",
            "The panel is deliberately small and directional; absent canonical sites cannot be substituted post hoc.",
            "Module residualisation adjusts each phosphosite for corresponding total protein before aggregation.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
