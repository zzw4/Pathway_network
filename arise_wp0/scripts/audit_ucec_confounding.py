#!/usr/bin/env python3
"""Audit whether UCEC PTEN/KRAS source labels collapse to known subtypes."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.stats import chi2_contingency, fisher_exact

from refine_functional_sources import API, matches


STUDY = "ucec_tcga_pan_can_atlas_2018"
GENES = {"PTEN": 5728, "KRAS": 3845}
CLASSES = ["PTEN_LOF", "KRAS_CANONICAL"]


def cramers_v(table: pd.DataFrame) -> float:
    if table.shape[0] < 2 or table.shape[1] < 2:
        return float("nan")
    chi2 = chi2_contingency(table, correction=False)[0]
    n = table.to_numpy().sum()
    return float(np.sqrt((chi2 / n) / min(table.shape[0] - 1, table.shape[1] - 1)))


def fetch_clinical(session: requests.Session, clinical_type: str) -> list[dict]:
    response = session.get(
        f"{API}/studies/{STUDY}/clinical-data",
        params={"clinicalDataType": clinical_type, "pageSize": 100000}, timeout=180,
    )
    response.raise_for_status()
    return response.json()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "ARISE-WP0/0.1"})

    samples = session.get(f"{API}/studies/{STUDY}/samples", params={"pageSize": 10000}, timeout=90).json()
    primary = [x for x in samples if "Primary" in x.get("sampleType", "")]
    sample_to_patient = {x["sampleId"]: x["patientId"] for x in primary}
    patients = sorted(set(sample_to_patient.values()))

    mutations = session.post(
        f"{API}/molecular-profiles/{STUDY}_mutations/mutations/fetch",
        params={"projection": "DETAILED"},
        json={"sampleIds": list(sample_to_patient), "entrezGeneIds": list(GENES.values())}, timeout=180,
    ).json()
    carriers = {name: set() for name in CLASSES}
    for record in mutations:
        sample = record.get("sampleId")
        if sample not in sample_to_patient:
            continue
        gene = (record.get("gene") or {}).get("hugoGeneSymbol", "")
        for name in CLASSES:
            if matches(name, gene, record):
                carriers[name].add(sample_to_patient[sample])

    sample_clin = defaultdict(dict)
    for row in fetch_clinical(session, "SAMPLE"):
        sid = row.get("sampleId")
        if sid in sample_to_patient:
            sample_clin[sample_to_patient[sid]][row["clinicalAttributeId"]] = row.get("value")
    patient_clin = defaultdict(dict)
    for row in fetch_clinical(session, "PATIENT"):
        patient_clin[row["patientId"]][row["clinicalAttributeId"]] = row.get("value")

    records = []
    for patient in patients:
        s = sample_clin[patient]
        p = patient_clin[patient]
        mantis = pd.to_numeric(s.get("MSI_SCORE_MANTIS"), errors="coerce")
        records.append({
            "patient": patient,
            "PTEN_LOF": int(patient in carriers["PTEN_LOF"]),
            "KRAS_CANONICAL": int(patient in carriers["KRAS_CANONICAL"]),
            "SUBTYPE": p.get("SUBTYPE"),
            "TUMOR_TYPE": s.get("TUMOR_TYPE"),
            "GRADE": s.get("GRADE"),
            "MSI_MANTIS": None if pd.isna(mantis) else float(mantis),
            "MSI_BINARY": None if pd.isna(mantis) else ("MSI" if mantis > 0.6 else "non-MSI"),
        })
    frame = pd.DataFrame(records)

    categorical = {}
    for source in CLASSES:
        categorical[source] = {}
        for covariate in ["SUBTYPE", "TUMOR_TYPE", "GRADE", "MSI_BINARY"]:
            use = frame.dropna(subset=[covariate])
            table = pd.crosstab(use[covariate], use[source])
            chi2, pvalue, _, _ = chi2_contingency(table) if table.shape[1] == 2 else (np.nan, np.nan, None, None)
            categorical[source][covariate] = {
                "n": int(len(use)), "table": {str(i): {str(k): int(v) for k, v in row.items()} for i, row in table.to_dict(orient="index").items()},
                "chi_square_p": float(pvalue), "cramers_v": cramers_v(table),
            }

    # Within-subtype 2x2 support: the core source-comparison design must not rely on subtype alone.
    within = {}
    for subtype, group in frame.dropna(subset=["SUBTYPE"]).groupby("SUBTYPE"):
        tab = pd.crosstab(group["PTEN_LOF"], group["KRAS_CANONICAL"]).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
        odds, pvalue = fisher_exact(tab.to_numpy())
        within[subtype] = {
            "n": int(len(group)), "cells": {"00": int(tab.loc[0, 0]), "01": int(tab.loc[0, 1]), "10": int(tab.loc[1, 0]), "11": int(tab.loc[1, 1])},
            "fisher_odds_ratio": None if not np.isfinite(odds) else float(odds), "fisher_p": float(pvalue),
        }

    payload = {
        "study": STUDY, "n_primary_patients": len(patients),
        "carrier_counts": {k: len(v) for k, v in carriers.items()},
        "covariate_completeness": {k: int(frame[k].notna().sum()) for k in ["SUBTYPE", "TUMOR_TYPE", "GRADE", "MSI_MANTIS"]},
        "source_covariate_associations": categorical,
        "within_subtype_pair_support": within,
        "subtype_counts": dict(Counter(frame["SUBTYPE"].dropna())),
        "interpretation_rule": "Strong covariate association is expected biologically; viability requires usable within-subtype contrasts or explicit out-of-subtype validation, not covariate independence.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
