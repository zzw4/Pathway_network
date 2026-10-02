#!/usr/bin/env python3
"""Validate a frozen KRAS-to-ERK phosphosite module across CPTAC cancers."""

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

from refine_functional_sources import API, matches


BIOLOGICAL_SITES = ["MAPK1_T185", "MAPK1_Y187", "MAPK3_T202", "MAPK3_Y204"]
GENES = {"KRAS": 3845, "MAPK1": 5594, "MAPK3": 5595}
COHORTS = {
    "UCEC": {
        "study": "ucec_cptac_2020", "phospho": "ucec_cptac_2020_phosphoproteome",
        "protein": "ucec_cptac_2020_protein_quantification",
        "category": "GENOMICS_SUBTYPE", "numeric": "PURITY_CANCER",
        "sites": {x: x for x in BIOLOGICAL_SITES},
    },
    "LUAD": {
        "study": "luad_cptac_2020", "phospho": "luad_cptac_2020_phosphoproteome",
        "protein": "luad_cptac_2020_protein_quantification",
        "category": "MRNA_EXPRESSION_SUBTYPE_TCGA", "numeric": "SMOKING_SCORE_WGS",
        "sites": {
            "MAPK1_T185_Y187_double": "NP_002736.3_2_2_185_187",
            "MAPK1_Y187": "NP_002736.3_1_1_187_187",
            "MAPK3_T202_Y204_double": "NP_002737.2_2_2_202_204",
            "MAPK3_Y204": "NP_002737.2_1_1_204_204",
        },
    },
    "COAD": {
        "study": "coad_cptac_2019", "phospho": "coad_cptac_2019_phosphoprotein_quantification",
        "protein": "coad_cptac_2019_protein_quantification",
        "category": "MSI_STATUS", "numeric": None,
        "sites": {x: x.replace("_", "_p", 1) for x in BIOLOGICAL_SITES},
    },
}


def clinical(session: requests.Session, study: str, kind: str) -> list[dict]:
    r = session.get(
        f"{API}/studies/{study}/clinical-data",
        params={"clinicalDataType": kind, "pageSize": 100000}, timeout=180,
    )
    r.raise_for_status()
    return r.json()


def directional_fit(frame: pd.DataFrame, outcome: str, protein: bool) -> dict:
    terms = ["KRAS_CANONICAL"]
    needed = [outcome, "KRAS_CANONICAL"]
    if protein:
        terms.append("TOTAL_PROTEIN")
        needed.append("TOTAL_PROTEIN")
    if frame["CATEGORY"].notna().mean() >= 0.8 and frame["CATEGORY"].nunique() > 1:
        counts = frame["CATEGORY"].value_counts()
        keep = set(counts[counts >= 3].index)
        frame = frame.copy()
        frame["CATEGORY"] = frame["CATEGORY"].where(frame["CATEGORY"].isin(keep), "OTHER")
        terms.append("C(CATEGORY)")
        needed.append("CATEGORY")
    if frame["NUMERIC"].notna().mean() >= 0.8:
        terms.append("NUMERIC")
        needed.append("NUMERIC")
    use = frame.dropna(subset=needed)
    if use["KRAS_CANONICAL"].sum() < 5:
        return {"n": len(use), "carriers": int(use["KRAS_CANONICAL"].sum()), "error": "insufficient carriers"}
    fit = smf.ols(f"{outcome} ~ {' + '.join(terms)}", data=use).fit(cov_type="HC3")
    beta, p2 = float(fit.params["KRAS_CANONICAL"]), float(fit.pvalues["KRAS_CANONICAL"])
    return {
        "n": int(len(use)), "carriers": int(use["KRAS_CANONICAL"].sum()), "beta": beta,
        "se": float(fit.bse["KRAS_CANONICAL"]), "ci95": [float(x) for x in fit.conf_int().loc["KRAS_CANONICAL"]],
        "p_two_sided": p2, "p_directional_positive": p2 / 2 if beta > 0 else 1 - p2 / 2,
        "covariates": terms,
    }


def audit(session: requests.Session, cancer: str, spec: dict) -> dict:
    study = spec["study"]
    samples = session.get(f"{API}/studies/{study}/samples", params={"pageSize": 10000}, timeout=90).json()
    sample_ids = [x["sampleId"] for x in samples]
    sample_to_patient = {x["sampleId"]: x["patientId"] for x in samples}
    mut = session.post(
        f"{API}/molecular-profiles/{study}_mutations/mutations/fetch", params={"projection": "DETAILED"},
        json={"sampleIds": sample_ids, "entrezGeneIds": [GENES["KRAS"]]}, timeout=180,
    )
    mut.raise_for_status()
    carriers = {x["patientId"] for x in mut.json() if matches("KRAS_CANONICAL", "KRAS", x)}

    sc, pc = defaultdict(dict), defaultdict(dict)
    for row in clinical(session, study, "SAMPLE"):
        sc[row["sampleId"]][row["clinicalAttributeId"]] = row.get("value")
    for row in clinical(session, study, "PATIENT"):
        pc[row["patientId"]][row["clinicalAttributeId"]] = row.get("value")

    phospho = {}
    for biological_site, feature_id in spec["sites"].items():
        r = session.get(f"{API}/generic-assay-data/{spec['phospho']}/generic-assay/{feature_id}", timeout=120)
        r.raise_for_status()
        phospho[biological_site] = {x["sampleId"]: pd.to_numeric(x.get("value"), errors="coerce") for x in r.json()}
    prot = session.post(
        f"{API}/molecular-profiles/{spec['protein']}/molecular-data/fetch", params={"projection": "DETAILED"},
        json={"sampleIds": sample_ids, "entrezGeneIds": [GENES["MAPK1"], GENES["MAPK3"]]}, timeout=180,
    )
    prot.raise_for_status()
    proteins = defaultdict(dict)
    for x in prot.json():
        gene = (x.get("gene") or {}).get("hugoGeneSymbol", "")
        proteins[gene][x["sampleId"]] = pd.to_numeric(x.get("value"), errors="coerce")

    rows = []
    for sid in sample_ids:
        patient = sample_to_patient[sid]
        merged = {**pc[patient], **sc[sid]}
        rows.append({
            "sample": sid, "KRAS_CANONICAL": int(patient in carriers),
            "CATEGORY": merged.get(spec["category"]),
            "NUMERIC": pd.to_numeric(merged.get(spec["numeric"]), errors="coerce") if spec["numeric"] else np.nan,
        })
    base = pd.DataFrame(rows)
    site_results, residual_cols = [], []
    for i, site in enumerate(spec["sites"]):
        gene = site.split("_")[0]
        frame = base.copy()
        frame["PHOSPHO"] = frame["sample"].map(phospho[site])
        frame["TOTAL_PROTEIN"] = frame["sample"].map(proteins[gene])
        result = directional_fit(frame, "PHOSPHO", protein=True)
        result["site"] = site
        result["feature_id"] = spec["sites"][site]
        site_results.append(result)
        avail = frame.dropna(subset=["PHOSPHO", "TOTAL_PROTEIN"])
        if len(avail) >= 30:
            fit = smf.ols("PHOSPHO ~ TOTAL_PROTEIN", data=avail).fit()
            resid = pd.Series(fit.resid.to_numpy(), index=avail["sample"].to_numpy())
            col = f"site_{i}"
            base[col] = base["sample"].map((resid - resid.mean()) / resid.std(ddof=0))
            residual_cols.append(col)
    valid = [i for i, x in enumerate(site_results) if "p_directional_positive" in x]
    if valid:
        q = multipletests([site_results[i]["p_directional_positive"] for i in valid], method="fdr_bh")[1]
        for i, value in zip(valid, q):
            site_results[i]["q_directional_within_cancer"] = float(value)
    base["ERK_MODULE"] = base[residual_cols].mean(axis=1, skipna=True)
    base.loc[base[residual_cols].notna().sum(axis=1) < 2, "ERK_MODULE"] = np.nan
    module = directional_fit(base, "ERK_MODULE", protein=False)
    passed = bool(module.get("beta", 0) > 0 and module.get("p_two_sided", 1) < 0.05 and any(x.get("q_directional_within_cancer", 1) < 0.1 for x in site_results))
    return {
        "cancer": cancer, "study": study, "samples": len(sample_ids), "kras_carriers": len(carriers),
        "category_covariate": spec["category"], "numeric_covariate": spec["numeric"],
        "site_results": site_results, "module_result": module, "passes_gate": passed,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "ARISE-WP0/0.1"})
    results = {cancer: audit(session, cancer, spec) for cancer, spec in COHORTS.items()}
    payload = {
        "frozen_biological_sites": BIOLOGICAL_SITES, "results": results,
        "gate": "Positive ERK module at two-sided P<0.05 plus at least one site at within-cancer directional FDR<0.10.",
        "overall_pass": sum(x["passes_gate"] for x in results.values()) >= 2,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
