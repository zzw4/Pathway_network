#!/usr/bin/env python3
"""CPTAC UCEC validation gates for CTNNB1-WNT and TP53-p53 routes."""

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


STUDY = "ucec_cptac_2020"
PHOSPHO = f"{STUDY}_phosphoproteome"
PROTEIN = f"{STUDY}_protein_quantification"
MRNA = f"{STUDY}_mrna"
GENES = {"CTNNB1": 1499, "TP53": 7157, "CDKN1A": 1026, "MDM2": 4193, "BAX": 581}
DESTRUCTION_SITES = ["CTNNB1_S33", "CTNNB1_S37", "CTNNB1_T41", "CTNNB1_S45"]


def ctnnb1_activating(record: dict) -> bool:
    if (record.get("gene") or {}).get("hugoGeneSymbol") != "CTNNB1":
        return False
    change = (record.get("proteinChange") or "").replace("p.", "").upper()
    return any(change.startswith(x) for x in ["D32", "S33", "G34", "S37", "T41", "S45"])


def clinical(session: requests.Session, kind: str) -> list[dict]:
    r = session.get(f"{API}/studies/{STUDY}/clinical-data", params={"clinicalDataType": kind, "pageSize": 100000}, timeout=180)
    r.raise_for_status()
    return r.json()


def fit(frame: pd.DataFrame, outcome: str, source: str, expected: int, total_protein: bool = False) -> dict:
    terms = [source]
    needed = [outcome, source]
    if total_protein:
        terms.append("TOTAL_PROTEIN")
        needed.append("TOTAL_PROTEIN")
    preliminary = frame.dropna(subset=needed)
    if preliminary["SUBTYPE"].notna().mean() >= 0.8 and preliminary["SUBTYPE"].nunique() > 1:
        terms.append("C(SUBTYPE)")
        needed.append("SUBTYPE")
    if preliminary["PURITY"].notna().mean() >= 0.8:
        terms.append("PURITY")
        needed.append("PURITY")
    use = frame.dropna(subset=needed)
    if len(use) < 20 or use[source].sum() < 3 or use[source].nunique() < 2:
        return {"n": len(use), "carriers": int(use[source].sum()) if len(use) else 0,
                "error": "insufficient evaluable groups", "expected_sign": expected}
    model = smf.ols(f"{outcome} ~ {' + '.join(terms)}", data=use).fit(cov_type="HC3")
    beta, p2 = float(model.params[source]), float(model.pvalues[source])
    return {
        "n": len(use), "carriers": int(use[source].sum()), "beta": beta,
        "ci95": [float(x) for x in model.conf_int().loc[source]], "p_two_sided": p2,
        "p_directional": p2 / 2 if beta * expected > 0 else 1 - p2 / 2,
        "expected_sign": expected,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "ARISE-WP0/0.1"})
    samples = session.get(f"{API}/studies/{STUDY}/samples", params={"pageSize": 10000}, timeout=90).json()
    ids = [x["sampleId"] for x in samples]
    sid_patient = {x["sampleId"]: x["patientId"] for x in samples}
    mut = session.post(
        f"{API}/molecular-profiles/{STUDY}_mutations/mutations/fetch", params={"projection": "DETAILED"},
        json={"sampleIds": ids, "entrezGeneIds": [GENES["CTNNB1"], GENES["TP53"]]}, timeout=180,
    )
    mut.raise_for_status()
    ctnnb1 = {x["patientId"] for x in mut.json() if ctnnb1_activating(x)}
    tp53 = {x["patientId"] for x in mut.json() if matches("TP53_DAMAGING", "TP53", x)}

    sc, pc = defaultdict(dict), defaultdict(dict)
    for row in clinical(session, "SAMPLE"):
        sc[row["sampleId"]][row["clinicalAttributeId"]] = row.get("value")
    for row in clinical(session, "PATIENT"):
        pc[row["patientId"]][row["clinicalAttributeId"]] = row.get("value")
    base = pd.DataFrame([{
        "sample": sid, "CTNNB1_ACTIVATING": int(sid_patient[sid] in ctnnb1), "TP53_DAMAGING": int(sid_patient[sid] in tp53),
        "SUBTYPE": sc[sid].get("GENOMICS_SUBTYPE") or pc[sid_patient[sid]].get("GENOMICS_SUBTYPE"),
        "PURITY": pd.to_numeric(sc[sid].get("PURITY_CANCER") or pc[sid_patient[sid]].get("PURITY_CANCER"), errors="coerce"),
    } for sid in ids])

    protein_response = session.post(
        f"{API}/molecular-profiles/{PROTEIN}/molecular-data/fetch", params={"projection": "DETAILED"},
        json={"sampleIds": ids, "entrezGeneIds": list(GENES.values())}, timeout=180,
    )
    protein_response.raise_for_status()
    proteins = defaultdict(dict)
    for x in protein_response.json():
        gene = (x.get("gene") or {}).get("hugoGeneSymbol", "")
        proteins[gene][x["sampleId"]] = pd.to_numeric(x.get("value"), errors="coerce")

    wnt_results = []
    frame = base.copy()
    frame["OUTCOME"] = frame["sample"].map(proteins["CTNNB1"])
    total_ctnnb1 = fit(frame, "OUTCOME", "CTNNB1_ACTIVATING", expected=1)
    total_ctnnb1["feature"] = "CTNNB1_total_protein"
    wnt_results.append(total_ctnnb1)
    residual_cols = []
    for i, site in enumerate(DESTRUCTION_SITES):
        r = session.get(f"{API}/generic-assay-data/{PHOSPHO}/generic-assay/{site}", timeout=120)
        r.raise_for_status()
        values = {x["sampleId"]: pd.to_numeric(x.get("value"), errors="coerce") for x in r.json()}
        frame = base.copy()
        frame["OUTCOME"] = frame["sample"].map(values)
        frame["TOTAL_PROTEIN"] = frame["sample"].map(proteins["CTNNB1"])
        result = fit(frame, "OUTCOME", "CTNNB1_ACTIVATING", expected=-1, total_protein=True)
        result["feature"] = site
        wnt_results.append(result)
        avail = frame.dropna(subset=["OUTCOME", "TOTAL_PROTEIN"])
        if len(avail) >= 30:
            m = smf.ols("OUTCOME ~ TOTAL_PROTEIN", data=avail).fit()
            residual = pd.Series(m.resid.to_numpy(), index=avail["sample"].to_numpy())
            col = f"wnt_{i}"
            base[col] = base["sample"].map((residual - residual.mean()) / residual.std(ddof=0))
            residual_cols.append(col)
    valid_wnt = [x for x in wnt_results if "p_directional" in x]
    q = multipletests([x["p_directional"] for x in valid_wnt], method="fdr_bh")[1]
    for row, value in zip(valid_wnt, q):
        row["q_directional"] = float(value)
    base["DESTRUCTION_PHOSPHO_MODULE"] = base[residual_cols].mean(axis=1, skipna=True)
    base.loc[base[residual_cols].notna().sum(axis=1) < 2, "DESTRUCTION_PHOSPHO_MODULE"] = np.nan
    wnt_module = fit(base, "DESTRUCTION_PHOSPHO_MODULE", "CTNNB1_ACTIVATING", expected=-1)

    p53_results = []
    for gene in ["CDKN1A", "MDM2", "BAX"]:
        frame = base.copy()
        frame["OUTCOME"] = frame["sample"].map(proteins[gene])
        result = fit(frame, "OUTCOME", "TP53_DAMAGING", expected=-1)
        result["feature"] = f"{gene}_total_protein"
        p53_results.append(result)
    valid_p53 = [x for x in p53_results if "p_directional" in x]
    q = multipletests([x["p_directional"] for x in valid_p53], method="fdr_bh")[1]
    for row, value in zip(valid_p53, q):
        row["q_directional"] = float(value)

    payload = {
        "study": STUDY, "samples": len(ids), "carrier_counts": {"CTNNB1_ACTIVATING": len(ctnnb1), "TP53_DAMAGING": len(tp53)},
        "wnt": {"features": wnt_results, "destruction_site_module": wnt_module,
                "passes": bool(total_ctnnb1.get("q_directional", 1) < 0.1 or (wnt_module.get("beta", 1) < 0 and wnt_module.get("p_two_sided", 1) < 0.05))},
        "p53": {"features": p53_results, "passes": any(x.get("q_directional", 1) < 0.1 for x in p53_results)},
        "notes": ["CTNNB1 activation predicts higher total beta-catenin and/or lower destruction-site phosphorylation.",
                  "TP53 conservative damaging class predicts lower abundance of canonical p53 target proteins."],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
