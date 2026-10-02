#!/usr/bin/env python3
"""Estimate source-specific overlap and effective sample size for ARISE WP0."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder


SOURCE_COLUMNS = {
    "PIK3CA_GOF": "pik3ca_gof",
    "TP53_LOF": "tp53_lof",
}


def ess(weights: np.ndarray) -> float:
    weights = np.asarray(weights, dtype=float)
    return float(weights.sum() ** 2 / np.square(weights).sum()) if weights.size and np.square(weights).sum() else 0.0


def nearest_opposite_hamming(x: np.ndarray, treatment: np.ndarray) -> dict:
    values = []
    for group in (0, 1):
        left = x[treatment == group]
        right = x[treatment != group]
        if len(left) == 0 or len(right) == 0:
            continue
        for row in left:
            values.append(float(np.abs(right - row).mean(axis=1).min()))
    arr = np.asarray(values)
    return {
        "mean": float(arr.mean()) if arr.size else None,
        "median": float(np.median(arr)) if arr.size else None,
        "fraction_exact_match": float(np.mean(arr == 0)) if arr.size else None,
        "fraction_distance_le_0_2": float(np.mean(arr <= 0.2)) if arr.size else None,
    }


def audit_source(df: pd.DataFrame, label: str, source_col: str) -> dict:
    other_source = next(col for col in SOURCE_COLUMNS.values() if col != source_col)
    covariates = ["er", "pr", "her2", other_source, "erbb2_amp", "ccnd1_amp"]
    work = df.dropna(subset=["er", "pr", "her2"]).copy()
    y = work[source_col].astype(int).to_numpy()
    x = work[covariates]

    categorical = ["er", "pr", "her2"]
    numeric = [other_source, "erbb2_amp", "ccnd1_amp"]
    model = make_pipeline(
        ColumnTransformer([
            ("receptor", OneHotEncoder(drop="if_binary", handle_unknown="ignore"), categorical),
            ("binary", "passthrough", numeric),
        ]),
        LogisticRegression(C=1.0, penalty="l2", max_iter=2000, class_weight=None),
    )
    folds = StratifiedKFold(n_splits=5, shuffle=True, random_state=20261001)
    propensity = cross_val_predict(model, x, y, cv=folds, method="predict_proba")[:, 1]
    clipped = np.clip(propensity, 0.02, 0.98)

    ate_weight = np.where(y == 1, 1.0 / clipped, 1.0 / (1.0 - clipped))
    overlap_weight = np.where(y == 1, 1.0 - clipped, clipped)
    transformed = pd.get_dummies(x, columns=categorical, dtype=int).to_numpy(dtype=float)

    treated_p = propensity[y == 1]
    control_p = propensity[y == 0]
    common_low = max(float(treated_p.min()), float(control_p.min()))
    common_high = min(float(treated_p.max()), float(control_p.max()))
    in_common = (propensity >= common_low) & (propensity <= common_high)

    return {
        "label": label,
        "n": int(len(y)),
        "carriers": int(y.sum()),
        "noncarriers": int((1 - y).sum()),
        "propensity_auc_oof": float(roc_auc_score(y, propensity)),
        "propensity_range_carriers": [float(treated_p.min()), float(treated_p.max())],
        "propensity_range_noncarriers": [float(control_p.min()), float(control_p.max())],
        "common_support_interval": [common_low, common_high],
        "fraction_in_common_support": float(in_common.mean()),
        "carriers_in_common_support": int(np.sum(in_common & (y == 1))),
        "noncarriers_in_common_support": int(np.sum(in_common & (y == 0))),
        "ate_ess_carriers": ess(ate_weight[y == 1]),
        "ate_ess_noncarriers": ess(ate_weight[y == 0]),
        "ate_max_weight": float(ate_weight.max()),
        "overlap_ess_carriers": ess(overlap_weight[y == 1]),
        "overlap_ess_noncarriers": ess(overlap_weight[y == 0]),
        "nearest_opposite_hamming": nearest_opposite_hamming(transformed, y),
        "covariates": covariates,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    df = pd.read_csv(args.manifest, sep="\t")
    df = df[df["strict_four_modality"] == 1].copy()
    result = {
        "analysis": "Exploratory WP0 overlap audit; receptor and co-alteration covariates only.",
        "sources": {
            label: audit_source(df, label, column) for label, column in SOURCE_COLUMNS.items()
        },
        "interaction_gate": {
            "required_min_per_cell": 20,
            "complete_receptor_cells": {
                f"pik3ca_{p}_tp53_{t}": int(((df["pik3ca_gof"] == p) & (df["tp53_lof"] == t) & df[["er", "pr", "her2"]].notna().all(axis=1)).sum())
                for p in (0, 1) for t in (0, 1)
            },
        },
        "limitations": [
            "This is a support diagnostic, not a causal-effect estimate.",
            "Purity, ploidy, histology, grade, age and technical covariates are not yet joined.",
            "Thresholds must be rerun on the frozen GDC-harmonized cohort.",
        ],
    }
    cells = result["interaction_gate"]["complete_receptor_cells"]
    result["interaction_gate"]["passes_20_per_cell"] = all(value >= 20 for value in cells.values())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
