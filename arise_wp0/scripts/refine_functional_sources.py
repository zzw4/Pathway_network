#!/usr/bin/env python3
"""Refine candidate alteration sources into conservative functional classes."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import requests


API = "https://www.cbioportal.org/api"
LOF = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Nonstop_Mutation"}
TP53_HOTSPOTS = {"R175H", "R248Q", "R248W", "R273H", "R273C", "R282W"}

CANCERS = {
    "UCEC": {
        "study": "ucec_tcga_pan_can_atlas_2018",
        "genes": {"PTEN": 5728, "PIK3CA": 5290, "TP53": 7157, "PIK3R1": 5295, "KRAS": 3845},
        "classes": ["PTEN_LOF", "PIK3CA_CANONICAL", "TP53_DAMAGING", "PIK3R1_LOF", "KRAS_CANONICAL"],
    },
    "LUAD": {
        "study": "luad_tcga_pan_can_atlas_2018",
        "genes": {"KRAS": 3845, "EGFR": 1956, "STK11": 6794, "KEAP1": 9817, "TP53": 7157},
        "classes": ["KRAS_CANONICAL", "EGFR_ACTIVATING", "STK11_LOF", "KEAP1_LOF", "TP53_DAMAGING"],
    },
    "COADREAD": {
        "study": "coadread_tcga_pan_can_atlas_2018",
        "genes": {"APC": 324, "KRAS": 3845, "BRAF": 673, "PIK3CA": 5290, "TP53": 7157},
        "classes": ["APC_LOF", "KRAS_CANONICAL", "BRAF_V600E", "PIK3CA_CANONICAL", "TP53_DAMAGING"],
    },
    "BRCA": {
        "study": "brca_tcga_pan_can_atlas_2018",
        "genes": {"PIK3CA": 5290, "TP53": 7157, "PTEN": 5728, "AKT1": 207},
        "classes": ["PIK3CA_CANONICAL", "TP53_DAMAGING", "PTEN_LOF", "AKT1_E17K"],
    },
}


def is_lof(record: dict) -> bool:
    return record.get("mutationType", "") in LOF


def protein(record: dict) -> str:
    return (record.get("proteinChange") or "").replace("p.", "").upper()


def matches(class_name: str, gene: str, record: dict) -> bool:
    change = protein(record)
    if class_name == "PTEN_LOF": return gene == "PTEN" and is_lof(record)
    if class_name == "PIK3R1_LOF": return gene == "PIK3R1" and is_lof(record)
    if class_name == "STK11_LOF": return gene == "STK11" and is_lof(record)
    if class_name == "KEAP1_LOF": return gene == "KEAP1" and is_lof(record)
    if class_name == "APC_LOF": return gene == "APC" and is_lof(record)
    if class_name == "PIK3CA_CANONICAL": return gene == "PIK3CA" and change in {"E542K", "E545K", "H1047R", "H1047L"}
    if class_name == "KRAS_CANONICAL": return gene == "KRAS" and (change.startswith("G12") or change.startswith("G13") or change.startswith("Q61"))
    if class_name == "BRAF_V600E": return gene == "BRAF" and change == "V600E"
    if class_name == "AKT1_E17K": return gene == "AKT1" and change == "E17K"
    if class_name == "EGFR_ACTIVATING":
        return gene == "EGFR" and (change == "L858R" or ("DEL" in change and any(x in change for x in ("E746", "L747", "S752"))))
    if class_name == "TP53_DAMAGING": return gene == "TP53" and (is_lof(record) or change in TP53_HOTSPOTS)
    return False


def screen(session: requests.Session, cancer: str, spec: dict) -> dict:
    study = spec["study"]
    samples = session.get(f"{API}/studies/{study}/samples", params={"pageSize": 10000}, timeout=90).json()
    primary = [x for x in samples if "Primary" in x.get("sampleType", "")]
    sample_to_patient = {x["sampleId"]: x["patientId"] for x in primary}
    patients = set(sample_to_patient.values())
    response = session.post(
        f"{API}/molecular-profiles/{study}_mutations/mutations/fetch",
        params={"projection": "DETAILED"},
        json={"sampleIds": list(sample_to_patient), "entrezGeneIds": list(spec["genes"].values())},
        timeout=180,
    )
    response.raise_for_status()
    records = response.json()
    carriers = {name: set() for name in spec["classes"]}
    for record in records:
        sample = record.get("sampleId")
        if sample not in sample_to_patient:
            continue
        gene = (record.get("gene") or {}).get("hugoGeneSymbol", "")
        for class_name in spec["classes"]:
            if matches(class_name, gene, record):
                carriers[class_name].add(sample_to_patient[sample])
    pairs = []
    for first, second in itertools.combinations(spec["classes"], 2):
        a, b = carriers[first], carriers[second]
        cells = {
            "00": len(patients - a - b), "10": len(a - b),
            "01": len(b - a), "11": len(a & b),
        }
        pairs.append({"first": first, "second": second, "cells": cells, "min_cell": min(cells.values()), "passes_20_each": min(cells.values()) >= 20})
    pairs.sort(key=lambda x: (x["passes_20_each"], x["min_cell"]), reverse=True)
    return {
        "cancer": cancer,
        "study": study,
        "primary_patients": len(patients),
        "functional_class_counts": {key: len(value) for key, value in carriers.items()},
        "pairs": pairs,
        "passing_pairs": sum(x["passes_20_each"] for x in pairs),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "ARISE-WP0/0.1"})
    result = {cancer: screen(session, cancer, spec) for cancer, spec in CANCERS.items()}
    payload = {
        "definitions": {
            "PIK3CA_CANONICAL": ["E542K", "E545K", "H1047R", "H1047L"],
            "KRAS_CANONICAL": ["G12*", "G13*", "Q61*"],
            "TP53_DAMAGING": sorted(LOF) + sorted(TP53_HOTSPOTS),
            "EGFR_ACTIVATING": ["L858R", "canonical exon-19 deletions"],
            "LOF": sorted(LOF),
        },
        "results": result,
        "limitations": [
            "Conservative functional definitions deliberately trade sensitivity for interpretability.",
            "PIK3R1 and KEAP1 missense loss-of-function are undercounted pending curated functional annotation.",
            "Counts are TCGA primary-patient mutation support, not CPTAC validation support.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({c: {"n": x["primary_patients"], "classes": x["functional_class_counts"], "passing_pairs": x["passing_pairs"], "best_pair": x["pairs"][0]} for c, x in result.items()}, indent=2))


if __name__ == "__main__":
    main()
