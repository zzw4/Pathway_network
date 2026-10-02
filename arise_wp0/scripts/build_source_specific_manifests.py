#!/usr/bin/env python3
"""Freeze source-specific signed paths and direct TF-target outputs for ARISE Gate A."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path


SOURCES = {
    "KRAS": {
        "alteration_definition": "canonical activating hotspot/GOF mutation",
        "injection_node": "KRAS",
        "injection_sign": 1,
        "locked_intermediates": [
            {"node": "MAPK1/MAPK3", "assay": "phosphoprotein", "readout": "ERK activation-site module"},
            {"node": "MAP2K1/MAP2K2", "assay": "phosphoprotein", "readout": "MEK activation-site module"},
        ],
    },
    "CTNNB1": {
        "alteration_definition": "activating exon-3 mutation at D32/S33/G34/S37/T41/S45",
        "injection_node": "CTNNB1",
        "injection_sign": 1,
        "locked_intermediates": [
            {"node": "CTNNB1", "assay": "protein", "readout": "total beta-catenin abundance"},
            {"node": "CTNNB1", "assay": "phosphoprotein", "readout": "beta-catenin phosphosite(s), if mapped"},
        ],
    },
}


def simple_paths(adjacency, source, max_edges):
    stack = [(source, [source], 1)]
    while stack:
        node, path, sign = stack.pop()
        if len(path) - 1 >= max_edges:
            continue
        for nxt, edge_sign in adjacency.get(node, []):
            if nxt in path:
                continue
            new_path = path + [nxt]
            new_sign = sign * edge_sign
            yield new_path, new_sign
            stack.append((nxt, new_path, new_sign))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", default="arise_wp0/outputs/wp4155_augmented_graph.json")
    parser.add_argument("--outputs", default="arise_wp0/outputs/wp4155_output_layer.json")
    parser.add_argument("--out", default="arise_wp0/outputs/source_specific_graph_manifest.json")
    parser.add_argument("--max-edges", type=int, default=8)
    args = parser.parse_args()

    graph = json.loads(Path(args.graph).read_text(encoding="utf-8"))
    output_layer = json.loads(Path(args.outputs).read_text(encoding="utf-8"))
    terminal_tfs = set(graph["audit"]["terminal_tfs"])

    adjacency = defaultdict(list)
    core_edges = []
    for edge in graph["core_edges"]:
        if edge.get("directed") and edge.get("sign") in (-1, 1):
            adjacency[edge["source"]].append((edge["target"], int(edge["sign"])))
            core_edges.append(edge)

    direct_edges = [
        e for e in output_layer["direct_tf_edges"]
        if e.get("sign") in (-1, 1) and e["source"] in terminal_tfs
    ]

    manifest = {
        "purpose": "ARISE Gate A: source-specific node-by-node signed propagation",
        "graph": graph.get("name"),
        "path_policy": {
            "simple_paths_only": True,
            "max_edges": args.max_edges,
            "neutral_edges_excluded": True,
            "conflicting_source_to_tf_signs": "exclude TF from primary directional target set",
            "conflicting_target_signs": "exclude target from primary directional target set",
            "decoder": "DoRothEA direct signed TF-target edges only; PROGENy excluded",
        },
        "sources": {},
    }

    for source, specification in SOURCES.items():
        paths_by_node = defaultdict(list)
        for path, path_sign in simple_paths(adjacency, source, args.max_edges):
            paths_by_node[path[-1]].append({"nodes": path, "sign": path_sign})

        reachable_nodes = sorted(paths_by_node)
        tf_paths = {tf: paths_by_node[tf] for tf in sorted(terminal_tfs & set(paths_by_node))}
        tf_signs = {tf: sorted({p["sign"] for p in paths}) for tf, paths in tf_paths.items()}
        unambiguous_tfs = {tf: signs[0] for tf, signs in tf_signs.items() if len(signs) == 1}
        ambiguous_tfs = {tf: signs for tf, signs in tf_signs.items() if len(signs) > 1}

        target_evidence = defaultdict(list)
        for edge in direct_edges:
            tf = edge["source"]
            if tf not in unambiguous_tfs:
                continue
            target_evidence[edge["target"]].append({
                "tf": tf,
                "source_to_tf_sign": unambiguous_tfs[tf],
                "tf_to_target_sign": int(edge["sign"]),
                "predicted_sign": unambiguous_tfs[tf] * int(edge["sign"]),
                "confidence": edge.get("confidence"),
                "resource": edge.get("resource"),
            })

        directional_targets = []
        conflicting_targets = []
        for target, evidence in sorted(target_evidence.items()):
            signs = sorted({e["predicted_sign"] for e in evidence})
            record = {"target": target, "predicted_signs": signs, "evidence": evidence}
            if len(signs) == 1:
                record["predicted_sign"] = signs[0]
                directional_targets.append(record)
            else:
                conflicting_targets.append(record)

        reachable_edge_pairs = {(p["nodes"][i], p["nodes"][i + 1]) for paths in paths_by_node.values() for p in paths for i in range(len(p["nodes"]) - 1)}
        reachable_edges = [e for e in core_edges if (e["source"], e["target"]) in reachable_edge_pairs]

        manifest["sources"][source] = {
            **specification,
            "reachable_nodes": reachable_nodes,
            "reachable_edges": reachable_edges,
            "tf_paths": tf_paths,
            "tf_path_signs": tf_signs,
            "unambiguous_tfs": unambiguous_tfs,
            "ambiguous_tfs": ambiguous_tfs,
            "directional_targets": directional_targets,
            "conflicting_targets": conflicting_targets,
            "audit": {
                "reachable_node_count": len(reachable_nodes),
                "reachable_edge_count": len(reachable_edges),
                "reachable_tf_count": len(tf_paths),
                "unambiguous_tf_count": len(unambiguous_tfs),
                "ambiguous_tf_count": len(ambiguous_tfs),
                "directional_target_count": len(directional_targets),
                "conflicting_target_count": len(conflicting_targets),
            },
        }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({source: data["audit"] for source, data in manifest["sources"].items()}, indent=2))


if __name__ == "__main__":
    main()
