#!/usr/bin/env python3
"""Screen the ten CPTAC cancer types for TCGA alteration-source support."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import requests


API = "https://www.cbioportal.org/api"

STUDIES = {
    "UCEC": ("ucec_tcga_pan_can_atlas_2018", {"PTEN": 5728, "PIK3CA": 5290, "PIK3R1": 5295, "TP53": 7157, "ARID1A": 8289, "KRAS": 3845}),
    "LUAD": ("luad_tcga_pan_can_atlas_2018", {"KRAS": 3845, "EGFR": 1956, "STK11": 6794, "KEAP1": 9817, "TP53": 7157, "NF1": 4763, "BRAF": 673}),
    "COADREAD": ("coadread_tcga_pan_can_atlas_2018", {"APC": 324, "KRAS": 3845, "TP53": 7157, "PIK3CA": 5290, "BRAF": 673, "SMAD4": 4089}),
    "BRCA": ("brca_tcga_pan_can_atlas_2018", {"PIK3CA": 5290, "TP53": 7157, "GATA3": 2625, "MAP3K1": 4214, "CDH1": 999, "PTEN": 5728, "AKT1": 207, "ERBB2": 2064}),
    "HNSC": ("hnsc_tcga_pan_can_atlas_2018", {"TP53": 7157, "PIK3CA": 5290, "FAT1": 2195, "NOTCH1": 4851, "CASP8": 841, "HRAS": 3265}),
    "GBM": ("gbm_tcga_pan_can_atlas_2018", {"EGFR": 1956, "PTEN": 5728, "NF1": 4763, "TP53": 7157, "IDH1": 3417, "PIK3CA": 5290}),
    "CCRCC": ("kirc_tcga_pan_can_atlas_2018", {"VHL": 7428, "PBRM1": 55193, "SETD2": 29072, "BAP1": 8314, "MTOR": 2475}),
    "LSCC": ("lusc_tcga_pan_can_atlas_2018", {"TP53": 7157, "NFE2L2": 4780, "KEAP1": 9817, "PIK3CA": 5290, "CDKN2A": 1029, "PTEN": 5728, "NOTCH1": 4851}),
    "OV": ("ov_tcga_pan_can_atlas_2018", {"TP53": 7157, "BRCA1": 672, "BRCA2": 675, "NF1": 4763, "RB1": 5925, "PIK3CA": 5290}),
    "PDAC": ("paad_tcga_pan_can_atlas_2018", {"KRAS": 3845, "TP53": 7157, "SMAD4": 4089, "CDKN2A": 1029, "RNF43": 54894, "ARID1A": 8289}),
}

SILENT = {"Silent", "Intron", "3'UTR", "5'UTR", "IGR", "RNA", "lincRNA"}
LOF = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Nonstop_Mutation"}


def get_json(session: requests.Session, url: str, **kwargs):
    response = session.get(url, timeout=90, **kwargs)
    response.raise_for_status()
    return response.json()


def post_json(session: requests.Session, url: str, payload: dict, **kwargs):
    response = session.post(url, json=payload, timeout=180, **kwargs)
    response.raise_for_status()
    return response.json()


def mutation_symbol(record: dict, id_to_symbol: dict[int, str]) -> str | None:
    gene = record.get("gene") or {}
    return gene.get("hugoGeneSymbol") or id_to_symbol.get(record.get("entrezGeneId"))


def screen_one(session: requests.Session, cancer: str, study: str, genes: dict[str, int]) -> dict:
    samples = get_json(session, f"{API}/studies/{study}/samples", params={"pageSize": 10000})
    primary = [x for x in samples if "Primary" in x.get("sampleType", "")]
    sample_ids = [x["sampleId"] for x in primary]
    sample_to_patient = {x["sampleId"]: x["patientId"] for x in primary}
    patients = set(sample_to_patient.values())

    profiles = get_json(session, f"{API}/studies/{study}/molecular-profiles")
    mutation_profile = next(x["molecularProfileId"] for x in profiles if x.get("molecularAlterationType") == "MUTATION_EXTENDED")
    mutations = post_json(
        session,
        f"{API}/molecular-profiles/{mutation_profile}/mutations/fetch",
        {"sampleIds": sample_ids, "entrezGeneIds": list(genes.values())},
        params={"projection": "DETAILED"},
    )
    id_to_symbol = {value: key for key, value in genes.items()}
    carriers_any = {gene: set() for gene in genes}
    carriers_lof = {gene: set() for gene in genes}
    carriers_hotspot = {gene: set() for gene in genes}
    for record in mutations:
        symbol = mutation_symbol(record, id_to_symbol)
        sample = record.get("sampleId")
        if symbol not in genes or sample not in sample_to_patient:
            continue
        patient = sample_to_patient[sample]
        variant = record.get("variantType") or record.get("mutationType") or record.get("variantClassification") or ""
        classification = record.get("mutationType") or record.get("variantClassification") or ""
        if variant not in SILENT and classification not in SILENT:
            carriers_any[symbol].add(patient)
        if classification in LOF or variant in LOF:
            carriers_lof[symbol].add(patient)
        if record.get("hotspot") is True:
            carriers_hotspot[symbol].add(patient)

    frequency = {}
    for gene in genes:
        frequency[gene] = {
            "any_nonsilent_n": len(carriers_any[gene]),
            "any_nonsilent_fraction": len(carriers_any[gene]) / len(patients) if patients else 0.0,
            "lof_n": len(carriers_lof[gene]),
            "hotspot_n": len(carriers_hotspot[gene]),
        }

    pairs = []
    for first, second in itertools.combinations(genes, 2):
        a, b = carriers_any[first], carriers_any[second]
        cells = {
            "00": len(patients - a - b),
            "10": len((patients & a) - b),
            "01": len((patients & b) - a),
            "11": len(patients & a & b),
        }
        pairs.append({
            "first": first,
            "second": second,
            "cells": cells,
            "min_cell": min(cells.values()),
            "passes_20_each": min(cells.values()) >= 20,
        })
    pairs.sort(key=lambda item: (item["passes_20_each"], item["min_cell"]), reverse=True)

    return {
        "cancer": cancer,
        "study": study,
        "primary_samples": len(primary),
        "primary_patients": len(patients),
        "candidate_frequencies": frequency,
        "top_pairs": pairs[:10],
        "number_pairs_passing_20_each": sum(pair["passes_20_each"] for pair in pairs),
        "notes": "Mutation support only; CNA, RNA intersection, covariate overlap and CPTAC carriers remain to be audited.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "ARISE-WP0/0.1"})
    results = {}
    for cancer, (study, genes) in STUDIES.items():
        print(f"Screening {cancer} ({study})", flush=True)
        results[cancer] = screen_one(session, cancer, study, genes)
    payload = {
        "source": "cBioPortal public REST API, TCGA PanCancer Atlas studies",
        "results": results,
        "limitations": [
            "This first pass uses candidate-gene nonsilent mutation support, not frozen functional classes.",
            "Passing mutation four-cell counts does not establish local covariate overlap or signalling separability.",
            "CPTAC carrier counts and non-RNA intermediate coverage are separate mandatory gates.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({k: {"n": v["primary_patients"], "passing_pairs": v["number_pairs_passing_20_each"], "best_pair": v["top_pairs"][0] if v["top_pairs"] else None} for k, v in results.items()}, indent=2))


if __name__ == "__main__":
    main()
