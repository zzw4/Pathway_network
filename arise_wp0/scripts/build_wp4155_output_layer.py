#!/usr/bin/env python3
"""Build a typed expression-output layer for the WP4155 signalling graph."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


MODULES = ["MAPK", "PI3K", "WNT", "p53", "EGFR"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", required=True, type=Path)
    parser.add_argument("--progeny", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--top", type=int, default=100)
    args = parser.parse_args()
    graph = json.loads(args.graph.read_text(encoding="utf-8"))

    # Direct TF regulation remains distinct from pathway-response footprints.
    direct_edges = []
    seen_direct = set()
    for edge in graph["tf_decoder"]:
        key = (edge["source"], edge["target"], int(edge["sign"]))
        if key in seen_direct:
            continue
        seen_direct.add(key)
        direct_edges.append({
            "source": edge["source"], "target": edge["target"], "sign": int(edge["sign"]),
            "weight": None, "evidence_type": "direct_tf_target", "resource": "DoRothEA",
            "confidence": edge["confidence"],
        })

    raw = pd.read_csv(args.progeny, sep="\t")
    model = raw.pivot_table(index=["genesymbol", "record_id"], columns="label", values="value", aggfunc="first").reset_index()
    model["p_value"] = pd.to_numeric(model["p_value"], errors="coerce")
    model["weight"] = pd.to_numeric(model["weight"], errors="coerce")
    footprint_edges = []
    module_targets = {}
    for module in MODULES:
        part = model[model["pathway"] == module].dropna(subset=["genesymbol", "weight", "p_value"])
        part = part.sort_values("p_value").drop_duplicates("genesymbol").head(args.top)
        module_targets[module] = part["genesymbol"].tolist()
        for row in part.itertuples():
            footprint_edges.append({
                "source": f"module:{module}", "target": row.genesymbol,
                "sign": 1 if row.weight > 0 else -1 if row.weight < 0 else 0,
                "weight": float(row.weight), "p_value": float(row.p_value),
                "evidence_type": "pathway_response_footprint", "resource": "PROGENy",
                "confidence": f"top_{args.top}",
            })

    by_target = defaultdict(list)
    for edge in direct_edges + footprint_edges:
        by_target[edge["target"]].append(edge)
    outputs = []
    for target, evidence in sorted(by_target.items()):
        types = sorted({x["evidence_type"] for x in evidence})
        signs = {x["sign"] for x in evidence if x["sign"] != 0}
        outputs.append({
            "gene": target, "evidence_types": types, "evidence_count": len(evidence),
            "has_direct_tf_evidence": "direct_tf_target" in types,
            "has_footprint_evidence": "pathway_response_footprint" in types,
            "sign_conflict_across_evidence": len(signs) > 1,
        })

    output_sets = {
        "anchor_both_evidence": sorted(x["gene"] for x in outputs if x["has_direct_tf_evidence"] and x["has_footprint_evidence"]),
        "direct_only": sorted(x["gene"] for x in outputs if x["has_direct_tf_evidence"] and not x["has_footprint_evidence"]),
        "footprint_only": sorted(x["gene"] for x in outputs if x["has_footprint_evidence"] and not x["has_direct_tf_evidence"]),
    }
    payload = {
        "graph_name": graph["name"], "module_nodes": [f"module:{x}" for x in MODULES],
        "direct_tf_edges": direct_edges, "footprint_edges": footprint_edges,
        "outputs": outputs, "output_sets": output_sets, "module_targets": module_targets,
        "audit": {
            "direct_tf_edges": len(direct_edges), "footprint_edges": len(footprint_edges),
            "unique_output_genes": len(outputs),
            "anchor_both_evidence": len(output_sets["anchor_both_evidence"]),
            "direct_only": len(output_sets["direct_only"]), "footprint_only": len(output_sets["footprint_only"]),
            "outputs_with_sign_conflict": sum(x["sign_conflict_across_evidence"] for x in outputs),
            "direct_edge_signs": dict(Counter(x["sign"] for x in direct_edges)),
            "footprint_edge_signs": dict(Counter(x["sign"] for x in footprint_edges)),
        },
        "usage_policy": {
            "primary_outputs": "anchor_both_evidence plus footprint_only genes for the relevant module",
            "mechanistic_claims": "Only direct_tf_target edges may be described as direct regulatory links.",
            "footprint_claims": "PROGENy edges indicate response association, not direct regulation.",
            "conflicts": "Conflicting signs remain explicit and are not silently averaged.",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(payload["audit"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
