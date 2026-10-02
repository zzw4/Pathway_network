#!/usr/bin/env python3
"""Cross-validated ridge feasibility baseline for WP4155 alteration-to-expression prediction."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeCV
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from audit_ucec_expression_signal import hgnc_map
from refine_functional_sources import API, matches


STUDY = "ucec_tcga_pan_can_atlas_2018"
SOURCE_GENES = {"KRAS": 3845, "CTNNB1": 1499, "TP53": 7157, "PIK3R1": 5295, "PTEN": 5728, "PIK3CA": 5290, "ERBB2": 2064, "FGFR2": 2263}


def ctnnb1_activating(record: dict) -> bool:
    change = (record.get("proteinChange") or "").replace("p.", "").upper()
    return any(change.startswith(x) for x in ["D32", "S33", "G34", "S37", "T41", "S45"])


def clinical(session: requests.Session, kind: str) -> list[dict]:
    r = session.get(f"{API}/studies/{STUDY}/clinical-data", params={"clinicalDataType": kind, "pageSize": 100000}, timeout=180)
    r.raise_for_status()
    return r.json()


def cross_validated_predictions(x: pd.DataFrame, y: np.ndarray, categorical: list[str], numeric: list[str], seed: int) -> np.ndarray:
    prediction = np.full_like(y, np.nan, dtype=float)
    folds = KFold(n_splits=5, shuffle=True, random_state=seed)
    for train, test in folds.split(x):
        prep = ColumnTransformer([
            ("cat", Pipeline([("impute", SimpleImputer(strategy="most_frequent")), ("onehot", OneHotEncoder(handle_unknown="ignore"))]), categorical),
            ("num", Pipeline([("impute", SimpleImputer(strategy="median")), ("scale", StandardScaler())]), numeric),
        ])
        model = Pipeline([("prep", prep), ("ridge", RidgeCV(alphas=np.logspace(-2, 4, 13)))])
        model.fit(x.iloc[train], y[train])
        prediction[test] = model.predict(x.iloc[test])
    return prediction


def summarise(y: np.ndarray, pred: np.ndarray, genes: list[str]) -> dict:
    per_gene = []
    for i, gene in enumerate(genes):
        keep = np.isfinite(y[:, i]) & np.isfinite(pred[:, i])
        value = r2_score(y[keep, i], pred[keep, i]) if keep.sum() >= 20 else np.nan
        corr = np.corrcoef(y[keep, i], pred[keep, i])[0, 1] if keep.sum() >= 20 else np.nan
        per_gene.append({"gene": gene, "r2": float(value), "pearson": float(corr)})
    values = np.array([x["r2"] for x in per_gene])
    return {"mean_r2": float(np.nanmean(values)), "median_r2": float(np.nanmedian(values)),
            "genes_positive_r2": int(np.nansum(values > 0)), "genes": len(per_gene), "per_gene": per_gene}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-layer", required=True, type=Path)
    parser.add_argument("--cache-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=20261001)
    parser.add_argument("--permutations", type=int, default=20)
    parser.add_argument("--dataset-output", type=Path)
    args = parser.parse_args()
    layer = json.loads(args.output_layer.read_text(encoding="utf-8"))
    output_symbols = sorted(set(layer["output_sets"]["anchor_both_evidence"] + layer["output_sets"]["footprint_only"]))
    mapping = hgnc_map(args.cache_dir / "hgnc_complete_set.txt")
    symbol_entrez = {x: mapping[x] for x in output_symbols if x in mapping}
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "ARISE-WP0/0.1"})
    samples = session.get(f"{API}/studies/{STUDY}/samples", params={"pageSize": 10000}, timeout=90).json()
    primary = [x for x in samples if "Primary" in x.get("sampleType", "")]
    sid_patient = {x["sampleId"]: x["patientId"] for x in primary}

    carriers = {x: set() for x in ["KRAS_CANONICAL", "CTNNB1_ACTIVATING", "TP53_DAMAGING", "PIK3R1_LOF", "PTEN_LOF", "PIK3CA_CANONICAL"]}
    mutation_response = session.post(
        f"{API}/molecular-profiles/{STUDY}_mutations/mutations/fetch", params={"projection": "DETAILED"},
        json={"sampleIds": list(sid_patient), "entrezGeneIds": list(SOURCE_GENES.values())}, timeout=240,
    )
    mutation_response.raise_for_status()
    for row in mutation_response.json():
        sid = row.get("sampleId")
        if sid not in sid_patient:
            continue
        gene = (row.get("gene") or {}).get("hugoGeneSymbol", "")
        patient = sid_patient[sid]
        for source in carriers:
            if source == "CTNNB1_ACTIVATING":
                if gene == "CTNNB1" and ctnnb1_activating(row): carriers[source].add(patient)
            elif matches(source, gene, row):
                carriers[source].add(patient)

    cna_response = session.post(
        f"{API}/molecular-profiles/{STUDY}_gistic/molecular-data/fetch",
        json={"sampleIds": list(sid_patient), "entrezGeneIds": list(SOURCE_GENES.values())}, timeout=240,
    )
    cna_response.raise_for_status()
    cna = defaultdict(dict)
    for row in cna_response.json():
        cna[row["sampleId"]][int(row["entrezGeneId"])] = pd.to_numeric(row.get("value"), errors="coerce")

    sc, pc = defaultdict(dict), defaultdict(dict)
    for row in clinical(session, "SAMPLE"):
        if row.get("sampleId") in sid_patient:
            sc[sid_patient[row["sampleId"]]][row["clinicalAttributeId"]] = row.get("value")
    for row in clinical(session, "PATIENT"):
        pc[row["patientId"]][row["clinicalAttributeId"]] = row.get("value")

    expr_response = session.post(
        f"{API}/molecular-profiles/{STUDY}_rna_seq_v2_mrna/molecular-data/fetch",
        json={"sampleIds": list(sid_patient), "entrezGeneIds": list(symbol_entrez.values())}, timeout=300,
    )
    expr_response.raise_for_status()
    expr = pd.DataFrame(expr_response.json()).pivot_table(index="sampleId", columns="entrezGeneId", values="value", aggfunc="first")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    genes = [g for g, e in symbol_entrez.items() if e in expr.columns]
    entrez_order = [symbol_entrez[g] for g in genes]
    expr = expr[entrez_order]
    expr.columns = genes
    expr = (expr - expr.mean()) / expr.std(ddof=0)

    rows = []
    for sid in expr.index:
        patient = sid_patient[sid]
        row = {"sample": sid, "SUBTYPE": pc[patient].get("SUBTYPE"), "TUMOR_TYPE": sc[patient].get("TUMOR_TYPE"),
               "MSI": pd.to_numeric(sc[patient].get("MSI_SCORE_MANTIS"), errors="coerce")}
        row.update({source: int(patient in members) for source, members in carriers.items()})
        for gene, entrez in SOURCE_GENES.items():
            row[f"CNA_{gene}"] = cna[sid].get(entrez, np.nan)
        rows.append(row)
    features = pd.DataFrame(rows).set_index("sample").loc[expr.index]
    keep_genes = expr.columns[expr.notna().mean() >= 0.9].tolist()
    y = expr[keep_genes].fillna(expr[keep_genes].mean()).to_numpy()

    categorical = ["SUBTYPE", "TUMOR_TYPE"]
    clinical_numeric = ["MSI"]
    source_numeric = list(carriers) + [f"CNA_{x}" for x in SOURCE_GENES]
    if args.dataset_output:
        clinical_frame = pd.concat([
            pd.get_dummies(features[categorical].fillna("MISSING").astype(str), prefix=categorical, dtype=float),
            features[clinical_numeric].apply(pd.to_numeric, errors="coerce").fillna(features[clinical_numeric].median()),
        ], axis=1)
        source_frame = features[source_numeric].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        args.dataset_output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            args.dataset_output,
            clinical=clinical_frame.to_numpy(dtype=np.float32), sources=source_frame.to_numpy(dtype=np.float32),
            y=y.astype(np.float32), sample_ids=np.asarray(features.index, dtype=str), genes=np.asarray(keep_genes, dtype=str),
            clinical_names=np.asarray(clinical_frame.columns, dtype=str), source_names=np.asarray(source_frame.columns, dtype=str),
        )
    pred_clinical = cross_validated_predictions(features, y, categorical, clinical_numeric, args.seed)
    pred_sources = cross_validated_predictions(features, y, categorical, clinical_numeric + source_numeric, args.seed)
    if args.dataset_output:
        np.savez_compressed(
            args.dataset_output,
            clinical=clinical_frame.to_numpy(dtype=np.float32), sources=source_frame.to_numpy(dtype=np.float32),
            y=y.astype(np.float32), clinical_oof=pred_clinical.astype(np.float32), source_ridge_oof=pred_sources.astype(np.float32),
            sample_ids=np.asarray(features.index, dtype=str), genes=np.asarray(keep_genes, dtype=str),
            clinical_names=np.asarray(clinical_frame.columns, dtype=str), source_names=np.asarray(source_frame.columns, dtype=str),
        )
    clinical_summary = summarise(y, pred_clinical, keep_genes)
    source_summary = summarise(y, pred_sources, keep_genes)
    deltas = []
    clinical_by_gene = {x["gene"]: x for x in clinical_summary["per_gene"]}
    for row in source_summary["per_gene"]:
        deltas.append({"gene": row["gene"], "delta_r2": row["r2"] - clinical_by_gene[row["gene"]]["r2"]})
    delta_values = np.array([x["delta_r2"] for x in deltas])
    gene_delta = {x["gene"]: x["delta_r2"] for x in deltas}
    module_increment = {}
    for module, targets in layer["module_targets"].items():
        values = np.array([gene_delta[x] for x in targets if x in gene_delta])
        module_increment[module] = {
            "genes": int(len(values)), "mean_delta_r2": float(values.mean()) if len(values) else None,
            "median_delta_r2": float(np.median(values)) if len(values) else None,
            "genes_delta_r2_gt_0.01": int((values > 0.01).sum()) if len(values) else 0,
        }

    # Shuffle all alteration/CNA columns jointly within molecular subtype, preserving their correlation structure.
    rng = np.random.default_rng(args.seed)
    null_mean_deltas = []
    groups = features["SUBTYPE"].fillna("MISSING").astype(str)
    for _ in range(args.permutations):
        shuffled = features.copy()
        for _, index_labels in groups.groupby(groups).groups.items():
            labels = list(index_labels)
            permuted = rng.permutation(labels)
            shuffled.loc[labels, source_numeric] = features.loc[permuted, source_numeric].to_numpy()
        null_pred = cross_validated_predictions(shuffled, y, categorical, clinical_numeric + source_numeric, args.seed)
        null_summary = summarise(y, null_pred, keep_genes)
        null_delta = np.mean([x["r2"] - clinical_by_gene[x["gene"]]["r2"] for x in null_summary["per_gene"]])
        null_mean_deltas.append(float(null_delta))
    empirical_p = (1 + sum(x >= float(delta_values.mean()) for x in null_mean_deltas)) / (1 + len(null_mean_deltas))
    payload = {
        "study": STUDY, "samples": len(features), "outputs": len(keep_genes),
        "carrier_counts": {x: len(y) for x, y in carriers.items()}, "source_features": source_numeric,
        "clinical_only": clinical_summary, "clinical_plus_alterations": source_summary,
        "increment": {"mean_delta_r2": float(delta_values.mean()), "median_delta_r2": float(np.median(delta_values)),
                      "genes_improved": int((delta_values > 0).sum()), "genes_delta_r2_gt_0.01": int((delta_values > 0.01).sum()),
                      "per_gene": deltas},
        "module_increment": module_increment,
        "subtype_stratified_permutation": {"permutations": args.permutations, "null_mean_deltas": null_mean_deltas,
                                           "null_mean": float(np.mean(null_mean_deltas)), "empirical_p": float(empirical_p)},
        "gate": "Alterations must improve mean held-out R2 and improve at least 10% of preregistered outputs by delta R2 > 0.01 before nonlinear modelling.",
    }
    payload["passes_gate"] = payload["increment"]["mean_delta_r2"] > 0 and payload["increment"]["genes_delta_r2_gt_0.01"] >= 0.1 * len(keep_genes) and empirical_p <= 0.05
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    compact = {k: v for k, v in payload.items() if k not in ["clinical_only", "clinical_plus_alterations", "increment"]}
    compact["clinical_only"] = {k: v for k, v in clinical_summary.items() if k != "per_gene"}
    compact["clinical_plus_alterations"] = {k: v for k, v in source_summary.items() if k != "per_gene"}
    compact["increment"] = {k: v for k, v in payload["increment"].items() if k != "per_gene"}
    print(json.dumps(compact, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
