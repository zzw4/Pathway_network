#!/usr/bin/env python3
"""Audit a cBioPortal-style TCGA-BRCA package using the Python standard library."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


REQUIRED = {
    "clinical": "data_clinical_sample.txt",
    "mutation": "data_mutations.txt",
    "cna": "data_cna.txt",
    "rna": "data_mrna_seq_v2_rsem.txt",
}

MISSING = {"", "na", "nan", "[not available]", "[not applicable]", "unknown", "indeterminate"}
PIK3CA_GOF = {"p.e542k", "p.e545k", "p.h1047r", "p.h1047l"}
TP53_LOF_CLASSES = {
    "nonsense_mutation",
    "frame_shift_del",
    "frame_shift_ins",
    "splice_site",
    "nonstop_mutation",
}


def patient_id(sample_id: str) -> str:
    value = sample_id.strip()
    return "-".join(value.split("-")[:3]) if value.upper().startswith("TCGA-") else value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normal_status(value: str) -> str | None:
    cleaned = value.strip().lower()
    if cleaned in MISSING:
        return None
    if cleaned == "positive":
        return "Positive"
    if cleaned == "negative":
        return "Negative"
    if cleaned == "equivocal":
        return "Equivocal"
    return None


def her2_status(ihc: str, fish: str) -> str | None:
    ihc_value = normal_status(ihc)
    fish_value = normal_status(fish)
    if ihc_value == "Positive":
        return "Positive"
    if ihc_value == "Negative":
        return "Negative"
    if fish_value in {"Positive", "Negative"}:
        return fish_value
    return None


def read_clinical(path: Path):
    rows = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        header = None
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if fields and fields[0] == "PATIENT_ID" and "SAMPLE_ID" in fields:
                header = fields
                break
        if header is None:
            raise ValueError(f"Could not find PATIENT_ID/SAMPLE_ID header in {path}")
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=header)
        for row in reader:
            if not row.get("PATIENT_ID") or not row.get("SAMPLE_ID"):
                continue
            if row.get("SAMPLE_TYPE", "").strip().lower() not in {"primary", "primary tumor", ""}:
                continue
            pid = patient_id(row["PATIENT_ID"])
            candidate = {
                "sample_id": row["SAMPLE_ID"],
                "er": normal_status(row.get("ER_STATUS_BY_IHC", "")),
                "pr": normal_status(row.get("PR_STATUS_BY_IHC", "")),
                "her2": her2_status(row.get("IHC_HER2", ""), row.get("HER2_FISH_STATUS", "")),
            }
            previous = rows.get(pid)
            if previous is None or sum(v is not None for k, v in candidate.items() if k != "sample_id") > sum(
                v is not None for k, v in previous.items() if k != "sample_id"
            ):
                rows[pid] = candidate
    return rows


def read_matrix_header(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    samples = [x for x in header[2:] if x]
    return set(samples), {patient_id(x) for x in samples}


def read_mutations(path: Path):
    all_patients = set()
    pik3ca = set()
    tp53_lof = set()
    tp53_missense = set()
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pid = patient_id(row.get("Tumor_Sample_Barcode", ""))
            if not pid:
                continue
            all_patients.add(pid)
            gene = row.get("Hugo_Symbol", "").upper()
            protein = row.get("HGVSp_Short", "").lower()
            klass = row.get("Variant_Classification", "").lower()
            if gene == "PIK3CA" and protein in PIK3CA_GOF:
                pik3ca.add(pid)
            if gene == "TP53" and klass in TP53_LOF_CLASSES:
                tp53_lof.add(pid)
            if gene == "TP53" and klass == "missense_mutation":
                tp53_missense.add(pid)
    return all_patients, pik3ca, tp53_lof, tp53_missense


def read_cna_sources(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        samples = header[2:]
        patients = {patient_id(x) for x in samples}
        amp = {"ERBB2": set(), "CCND1": set()}
        for row in reader:
            if not row:
                continue
            gene = row[0].upper()
            if gene not in amp:
                continue
            for sample, value in zip(samples, row[2:]):
                try:
                    if float(value) >= 2:
                        amp[gene].add(patient_id(sample))
                except ValueError:
                    pass
    return set(samples), patients, amp


def two_by_two(universe, first, second):
    return {
        "first0_second0": len(universe - first - second),
        "first1_second0": len((universe & first) - second),
        "first0_second1": len((universe & second) - first),
        "first1_second1": len(universe & first & second),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    paths = {key: args.study_dir / name for key, name in REQUIRED.items()}
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required files: " + ", ".join(missing))

    clinical = read_clinical(paths["clinical"])
    mutation_patients, pik3ca, tp53_lof, tp53_missense = read_mutations(paths["mutation"])
    _, cna_patients, amps = read_cna_sources(paths["cna"])
    _, rna_patients = read_matrix_header(paths["rna"])

    clinical_patients = set(clinical)
    strict = clinical_patients & mutation_patients & cna_patients & rna_patients
    receptor_complete = {
        pid for pid in strict if all(clinical[pid][key] is not None for key in ("er", "pr", "her2"))
    }

    result = {
        "source_study_dir": str(args.study_dir.resolve()),
        "counts": {
            "clinical_primary_patients": len(clinical_patients),
            "mutation_patients": len(mutation_patients),
            "cna_patients": len(cna_patients),
            "rna_patients": len(rna_patients),
            "strict_four_modality_patients": len(strict),
            "strict_with_complete_er_pr_her2": len(receptor_complete),
        },
        "strict_receptor_availability": {
            key: Counter(clinical[pid][key] or "Missing" for pid in strict)
            for key in ("er", "pr", "her2")
        },
        "sources_in_strict": {
            "PIK3CA_GOF_hotspot": len(strict & pik3ca),
            "TP53_LOF": len(strict & tp53_lof),
            "TP53_missense_separate": len(strict & tp53_missense),
            "ERBB2_high_level_amp": len(strict & amps["ERBB2"]),
            "CCND1_high_level_amp": len(strict & amps["CCND1"]),
        },
        "pik3ca_gof_by_tp53_lof_in_strict": two_by_two(strict, pik3ca, tp53_lof),
        "pik3ca_gof_by_tp53_lof_complete_receptors": two_by_two(receptor_complete, pik3ca, tp53_lof),
        "file_sha256": {key: sha256(path) for key, path in paths.items()},
        "limitations": [
            "Legacy cBioPortal/TCGA 2015-style package; not yet the frozen GDC-harmonized training set.",
            "Counts establish marginal support, not propensity overlap or local effective sample size.",
            "TP53 missense variants are intentionally not pooled with truncating/splice LOF.",
            "Clinical HER2 uses IHC with FISH resolution where available; unresolved equivocal values are missing.",
        ],
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "tcga_legacy_audit.json"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")

    tsv_path = args.output_dir / "patient_manifest.tsv"
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "patient_id", "strict_four_modality", "er", "pr", "her2", "pik3ca_gof",
            "tp53_lof", "tp53_missense", "erbb2_amp", "ccnd1_amp"
        ])
        for pid in sorted(clinical_patients | mutation_patients | cna_patients | rna_patients):
            status = clinical.get(pid, {})
            writer.writerow([
                pid, int(pid in strict), status.get("er", ""), status.get("pr", ""), status.get("her2", ""),
                int(pid in pik3ca), int(pid in tp53_lof), int(pid in tp53_missense),
                int(pid in amps["ERBB2"]), int(pid in amps["CCND1"]),
            ])

    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
