#!/usr/bin/env python3
"""Test whether UCEC PTEN/KRAS labels retain expected PROGENy RNA signals after covariate adjustment."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import statsmodels.formula.api as smf

from refine_functional_sources import API, matches


STUDY = "ucec_tcga_pan_can_atlas_2018"
PROFILE = f"{STUDY}_rna_seq_v2_mrna"
GENES = {"PTEN": 5728, "KRAS": 3845, "PIK3CA": 5290, "PIK3R1": 5295}


def clinical(session: requests.Session, kind: str) -> list[dict]:
    r = session.get(f"{API}/studies/{STUDY}/clinical-data", params={"clinicalDataType": kind, "pageSize": 100000}, timeout=180)
    r.raise_for_status()
    return r.json()


def get_progeny_model(cache: Path) -> pd.DataFrame:
    if not cache.exists():
        url = "https://omnipathdb.org/annotations?databases=PROGENy&format=tsv"
        response = requests.get(url, timeout=240)
        response.raise_for_status()
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_bytes(response.content)
    raw = pd.read_csv(cache, sep="\t")
    wide = raw.pivot_table(index=["genesymbol", "record_id"], columns="label", values="value", aggfunc="first").reset_index()
    wide["p_value"] = pd.to_numeric(wide["p_value"], errors="coerce")
    wide["weight"] = pd.to_numeric(wide["weight"], errors="coerce")
    return wide[wide["pathway"].isin(["PI3K", "MAPK"])].dropna(subset=["p_value", "weight"])


def hgnc_map(cache: Path) -> dict[str, int]:
    if not cache.exists():
        url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
        response = requests.get(url, timeout=180)
        response.raise_for_status()
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_bytes(response.content)
    data = pd.read_csv(cache, sep="\t", low_memory=False)
    data = data.dropna(subset=["symbol", "entrez_id"])
    return {str(r.symbol): int(r.entrez_id) for r in data.itertuples()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cache-dir", required=True, type=Path)
    parser.add_argument("--top", type=int, default=100)
    args = parser.parse_args()
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "ARISE-WP0/0.1"})

    samples = session.get(f"{API}/studies/{STUDY}/samples", params={"pageSize": 10000}, timeout=90).json()
    primary = [x for x in samples if "Primary" in x.get("sampleType", "")]
    sample_to_patient = {x["sampleId"]: x["patientId"] for x in primary}

    mutations = session.post(
        f"{API}/molecular-profiles/{STUDY}_mutations/mutations/fetch", params={"projection": "DETAILED"},
        json={"sampleIds": list(sample_to_patient), "entrezGeneIds": list(GENES.values())}, timeout=180,
    ).json()
    carriers = {"PTEN_LOF": set(), "KRAS_CANONICAL": set(), "PIK3CA_CANONICAL": set(), "PIK3R1_LOF": set()}
    for rec in mutations:
        sid = rec.get("sampleId")
        gene = (rec.get("gene") or {}).get("hugoGeneSymbol", "")
        if sid not in sample_to_patient:
            continue
        for source in carriers:
            if matches(source, gene, rec):
                carriers[source].add(sample_to_patient[sid])

    cna_response = session.post(
        f"{API}/molecular-profiles/{STUDY}_gistic/molecular-data/fetch",
        json={"sampleIds": list(sample_to_patient), "entrezGeneIds": [GENES["PTEN"]]}, timeout=180,
    )
    cna_response.raise_for_status()
    pten_cna = {row["sampleId"]: int(float(row["value"])) for row in cna_response.json() if row.get("value") is not None}
    pten_deepdel = {sample_to_patient[sid] for sid, value in pten_cna.items() if sid in sample_to_patient and value == -2}
    pten_anydel = {sample_to_patient[sid] for sid, value in pten_cna.items() if sid in sample_to_patient and value <= -1}
    carriers["PTEN_DEEPDEL"] = pten_deepdel
    carriers["PTEN_FUNCTIONAL"] = carriers["PTEN_LOF"] | pten_deepdel
    carriers["PTEN_BIALLELIC_PROXY"] = pten_deepdel | (carriers["PTEN_LOF"] & pten_anydel)

    sc, pc = defaultdict(dict), defaultdict(dict)
    for row in clinical(session, "SAMPLE"):
        if row.get("sampleId") in sample_to_patient:
            sc[sample_to_patient[row["sampleId"]]][row["clinicalAttributeId"]] = row.get("value")
    for row in clinical(session, "PATIENT"):
        pc[row["patientId"]][row["clinicalAttributeId"]] = row.get("value")

    model = get_progeny_model(args.cache_dir / "progeny_omnipath.tsv")
    mapping = hgnc_map(args.cache_dir / "hgnc_complete_set.txt")
    selected = []
    for pathway, group in model.groupby("pathway"):
        group = group.sort_values("p_value").drop_duplicates("genesymbol").head(args.top).copy()
        group["entrez"] = group["genesymbol"].map(mapping)
        selected.append(group.dropna(subset=["entrez"]))
    selected = pd.concat(selected, ignore_index=True)
    selected["entrez"] = selected["entrez"].astype(int)

    response = session.post(
        f"{API}/molecular-profiles/{PROFILE}/molecular-data/fetch",
        json={"sampleIds": list(sample_to_patient), "entrezGeneIds": sorted(selected["entrez"].unique().tolist())}, timeout=240,
    )
    response.raise_for_status()
    expr = pd.DataFrame(response.json())
    matrix = expr.pivot_table(index="sampleId", columns="entrezGeneId", values="value", aggfunc="first")
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    z = (matrix - matrix.mean()) / matrix.std(ddof=0)

    scores = {}
    used = {}
    for pathway, group in selected.groupby("pathway"):
        group = group[group["entrez"].isin(z.columns)].drop_duplicates("entrez")
        weights = group.set_index("entrez")["weight"].astype(float)
        block = z[weights.index]
        scores[pathway] = block.mul(weights, axis=1).sum(axis=1, min_count=max(20, len(weights)//2)) / np.sqrt(np.square(weights).sum())
        used[pathway] = {"genes": int(len(weights)), "expression_samples": int(scores[pathway].notna().sum())}

    score_frame = pd.DataFrame(scores)
    rows = []
    for sid, values in score_frame.iterrows():
        patient = sample_to_patient.get(sid)
        if patient is None:
            continue
        s, p = sc[patient], pc[patient]
        rows.append({
            "sample": sid, "patient": patient, **values.to_dict(),
            **{source: int(patient in members) for source, members in carriers.items()},
            "SUBTYPE": p.get("SUBTYPE"), "TUMOR_TYPE": s.get("TUMOR_TYPE"), "GRADE": s.get("GRADE"),
            "MSI_MANTIS": pd.to_numeric(s.get("MSI_SCORE_MANTIS"), errors="coerce"),
        })
    frame = pd.DataFrame(rows).dropna(subset=["SUBTYPE", "TUMOR_TYPE", "GRADE", "MSI_MANTIS"])

    results = {}
    formula_covars = "C(SUBTYPE) + C(TUMOR_TYPE) + C(GRADE) + MSI_MANTIS"
    source_names = ["PTEN_LOF", "PTEN_DEEPDEL", "PTEN_FUNCTIONAL", "PTEN_BIALLELIC_PROXY", "PIK3CA_CANONICAL", "PIK3R1_LOF", "KRAS_CANONICAL"]
    for pathway in ["PI3K", "MAPK"]:
        results[pathway] = {}
        for source in source_names:
            unadjusted = smf.ols(f"{pathway} ~ {source}", data=frame).fit(cov_type="HC3")
            other = "KRAS_CANONICAL" if source != "KRAS_CANONICAL" else "PTEN_FUNCTIONAL"
            adjusted = smf.ols(f"{pathway} ~ {source} + {other} + {formula_covars}", data=frame).fit(cov_type="HC3")
            results[pathway][source] = {
                "unadjusted_beta": float(unadjusted.params[source]), "unadjusted_p": float(unadjusted.pvalues[source]),
                "adjusted_beta": float(adjusted.params[source]), "adjusted_p": float(adjusted.pvalues[source]),
                "adjusted_ci95": [float(x) for x in adjusted.conf_int().loc[source].tolist()],
            }

    payload = {
        "study": STUDY, "n_complete_cases": int(len(frame)), "top_genes_requested_per_pathway": args.top,
        "carrier_counts": {source: int(frame[source].sum()) for source in source_names},
        "footprint_coverage": used, "results": results,
        "expected_checks": {"PTEN loss definitions": "positive PI3K activity", "KRAS_CANONICAL": "positive MAPK activity"},
        "warning": "PROGENy is an RNA footprint feasibility check, not independent validation and not the proposed model endpoint.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
