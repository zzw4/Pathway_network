#!/usr/bin/env python3
"""Audit candidate pathway intermediate coverage in CPTAC profiles."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import requests


API = "https://www.cbioportal.org/api"

GENES = {
    "AKT1": 207, "AKT2": 208, "MTOR": 2475, "RPS6KB1": 6198, "RPS6": 6194,
    "FOXO1": 2308, "FOXO3": 2309, "GSK3B": 2932, "TSC2": 7249,
    "EIF4EBP1": 1978, "AKT1S1": 84335, "RAF1": 5894, "BRAF": 673,
    "MAP2K1": 5604, "MAP2K2": 5605, "MAPK1": 5594, "MAPK3": 5595,
    "RPS6KA1": 6195, "ELK1": 2002, "CTNNB1": 1499, "AXIN1": 8312,
    "TCF7L2": 6934, "DVL1": 1855, "DVL2": 1856, "EGFR": 1956,
    "ERBB2": 2064, "SHC1": 6464, "GRB2": 2885, "SOS1": 6654, "GAB1": 2549,
    "TP53": 7157, "MDM2": 4193, "CDKN1A": 1026, "BAX": 581, "BBC3": 27113,
}

COHORTS = {
    "UCEC": {
        "study": "ucec_cptac_2020",
        "phospho": "ucec_cptac_2020_phosphoproteome",
        "protein": "ucec_cptac_2020_protein_quantification",
        "sources": {
            "PTEN_LOF": ["AKT1", "AKT2", "MTOR", "RPS6KB1", "RPS6", "FOXO1", "FOXO3", "GSK3B", "TSC2", "EIF4EBP1", "AKT1S1"],
            "KRAS_CANONICAL": ["RAF1", "BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "RPS6KA1", "ELK1"],
        },
    },
    "COADREAD": {
        "study": "coad_cptac_2019",
        "phospho": "coad_cptac_2019_phosphoprotein_quantification",
        "protein": "coad_cptac_2019_protein_quantification",
        "sources": {
            "APC_LOF": ["CTNNB1", "GSK3B", "AXIN1", "TCF7L2", "DVL1", "DVL2"],
            "KRAS_CANONICAL": ["RAF1", "BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "RPS6KA1", "ELK1"],
        },
    },
    "LUAD": {
        "study": "luad_cptac_2020",
        "phospho": "luad_cptac_2020_phosphoproteome",
        "protein": "luad_cptac_2020_protein_quantification",
        "sources": {
            "KRAS_CANONICAL": ["RAF1", "BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "RPS6KA1", "ELK1"],
            "EGFR_ACTIVATING": ["EGFR", "ERBB2", "SHC1", "GRB2", "SOS1", "GAB1", "AKT1", "MAPK1", "MAPK3"],
        },
    },
    "BRCA": {
        "study": "brca_cptac_2020",
        "phospho": "brca_cptac_2020_phosphoproteome",
        "protein": "brca_cptac_2020_protein_quantification",
        "sources": {
            "PIK3CA_CANONICAL": ["AKT1", "AKT2", "MTOR", "RPS6KB1", "RPS6", "FOXO1", "FOXO3", "GSK3B", "TSC2", "EIF4EBP1", "AKT1S1"],
            "TP53_DAMAGING": ["TP53", "MDM2", "CDKN1A", "BAX", "BBC3"],
        },
    },
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "ARISE-WP0/0.1"})
    results = {}
    for cancer, spec in COHORTS.items():
        samples = session.get(f"{API}/studies/{spec['study']}/samples", params={"pageSize": 10000}, timeout=90).json()
        sample_ids = [x["sampleId"] for x in samples]
        meta_response = session.get(
            f"{API}/generic-assay-meta/{spec['phospho']}",
            params={"projection": "DETAILED", "pageSize": 100000}, timeout=180,
        )
        meta_response.raise_for_status()
        phospho_symbols = Counter()
        for item in meta_response.json():
            props = item.get("genericEntityMetaProperties") or {}
            symbol = props.get("GENE_SYMBOL")
            if symbol:
                phospho_symbols[symbol.upper()] += 1

        all_symbols = sorted({gene for genes in spec["sources"].values() for gene in genes})
        protein_response = session.post(
            f"{API}/molecular-profiles/{spec['protein']}/molecular-data/fetch",
            params={"projection": "DETAILED"},
            json={"sampleIds": sample_ids, "entrezGeneIds": [GENES[g] for g in all_symbols]},
            timeout=180,
        )
        protein_response.raise_for_status()
        protein_samples = {gene: set() for gene in all_symbols}
        for item in protein_response.json():
            gene = (item.get("gene") or {}).get("hugoGeneSymbol", "").upper()
            value = item.get("value")
            if gene in protein_samples and value is not None:
                protein_samples[gene].add(item.get("sampleId"))

        source_result = {}
        for source, genes in spec["sources"].items():
            source_result[source] = {
                "candidate_intermediates": len(genes),
                "genes_with_phosphosites": sum(phospho_symbols[g] > 0 for g in genes),
                "total_phosphosite_features": sum(phospho_symbols[g] for g in genes),
                "phosphosite_features_by_gene": {g: phospho_symbols[g] for g in genes},
                "genes_with_protein_in_at_least_80pct_samples": sum(len(protein_samples[g]) >= 0.8 * len(sample_ids) for g in genes),
                "protein_nonmissing_samples_by_gene": {g: len(protein_samples[g]) for g in genes},
            }
        results[cancer] = {
            "study": spec["study"], "samples": len(sample_ids), "sources": source_result,
        }
        print(cancer, json.dumps(source_result), flush=True)

    payload = {
        "results": results,
        "limitations": [
            "Coverage is not evidence of source association or correct phosphosite direction.",
            "Candidate gene lists are a first mechanistic audit and must be frozen from curated signed resources.",
            "Site localisation and total-protein-adjusted usability are not yet tested.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    main()
