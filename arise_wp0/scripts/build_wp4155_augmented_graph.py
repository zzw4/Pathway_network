#!/usr/bin/env python3
"""Build a complete WP4155 scaffold augmented only by curated within-node-set edges."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict, deque
from pathlib import Path

import pandas as pd


ALTERATION_SOURCES = ["KRAS", "CTNNB1", "TP53", "PIK3R1", "PTEN", "PIK3CA", "ERBB2", "FGFR2"]
CORE_BRIDGES = {
    ("ERBB2", "GRB2"), ("ERBB2", "PIK3CA"), ("ERBB2", "PIK3R1"),
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--omnipath", required=True, type=Path)
    parser.add_argument("--dorothea", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))

    node_by_id = {x["id"]: x for x in manifest["nodes"]}
    group_by_id = {x["id"]: x for x in manifest["groups"]}
    genes = sorted({x["label"] for x in manifest["nodes"] if x["type"] == "GeneProduct"})
    gene_set = set(genes)
    signalling_entities = gene_set | {"PIP3"}

    def resolve_genes(entity_id: str | None) -> list[str]:
        if not entity_id:
            return []
        if entity_id in node_by_id:
            node = node_by_id[entity_id]
            return [node["label"]] if node["label"] in signalling_entities else []
        if entity_id in group_by_id:
            return sorted({node_by_id[x]["label"] for x in group_by_id[entity_id]["members"] if node_by_id[x]["label"] in signalling_entities})
        return []

    edge_map = {}
    for interaction in manifest["interactions"]:
        sources, targets = resolve_genes(interaction["source"]), resolve_genes(interaction["target"])
        for source in sources:
            for target in targets:
                if source == target:
                    continue
                key = (source, target, interaction["sign"])
                edge_map[key] = {
                    "source": source, "target": target, "sign": interaction["sign"],
                    "directed": interaction["directed"], "provenance": ["WikiPathways:WP4155"],
                    "wp_arrowhead": interaction["arrowhead"],
                }

    omni = pd.read_csv(args.omnipath, sep="\t", low_memory=False)
    omni = omni[(omni["is_directed"] == True) & (omni["consensus_direction"] == True)].copy()  # noqa: E712
    omni = omni[omni["source_genesymbol"].isin(gene_set) & omni["target_genesymbol"].isin(gene_set)]
    omni = omni[(omni["consensus_stimulation"] == True) | (omni["consensus_inhibition"] == True)]  # noqa: E712
    for row in omni.itertuples():
        signs = []
        if row.consensus_stimulation:
            signs.append(1)
        if row.consensus_inhibition:
            signs.append(-1)
        for sign in signs:
            key = (row.source_genesymbol, row.target_genesymbol, sign)
            provenance = [f"OmniPath:{x}" for x in str(row.sources).split(";") if x]
            if key in edge_map:
                edge_map[key]["provenance"] = sorted(set(edge_map[key]["provenance"] + provenance))
            else:
                edge_map[key] = {
                    "source": row.source_genesymbol, "target": row.target_genesymbol, "sign": sign,
                    "directed": True, "provenance": sorted(set(provenance)), "wp_arrowhead": None,
                }

    edges = list(edge_map.values())
    core_edges = [x for x in edges if "WikiPathways:WP4155" in x["provenance"] or (x["source"], x["target"]) in CORE_BRIDGES]
    ambiguous_pairs = defaultdict(set)
    for edge in edges:
        ambiguous_pairs[(edge["source"], edge["target"])].add(edge["sign"])

    regulons = pd.read_csv(args.dorothea, sep="\t", low_memory=False)
    regulons = regulons[regulons["dorothea_level"].astype(str).str.contains("A|B", regex=True)].copy()
    terminal_tfs = sorted(gene_set & set(regulons["source_genesymbol"]))
    decoder = []
    for row in regulons[regulons["source_genesymbol"].isin(terminal_tfs)].itertuples():
        sign = -1 if bool(row.consensus_inhibition) and not bool(row.consensus_stimulation) else 1 if bool(row.consensus_stimulation) and not bool(row.consensus_inhibition) else 0
        decoder.append({"source": row.source_genesymbol, "target": row.target_genesymbol, "sign": sign,
                        "confidence": row.dorothea_level, "provenance": "DoRothEA"})

    adjacency = defaultdict(set)
    for edge in core_edges:
        if edge["directed"]:
            adjacency[edge["source"]].add(edge["target"])

    reachability = {}
    for source in ALTERATION_SOURCES:
        distance = {source: 0}
        queue = deque([source])
        while queue:
            current = queue.popleft()
            for nxt in adjacency.get(current, set()):
                if nxt not in distance:
                    distance[nxt] = distance[current] + 1
                    queue.append(nxt)
        reachable_tfs = sorted(set(distance) & set(terminal_tfs))
        reachability[source] = {
            "reachable_nodes": len(distance), "reachable_tfs": reachable_tfs,
            "max_distance": max(distance.values()) if distance else None,
        }

    payload = {
        "name": "WP4155 complete augmented signalling graph", "wp_version": manifest["pathway"]["version"],
        "nodes": ([{"id": x, "type": "gene"} for x in genes] + [{"id": "PIP3", "type": "metabolite"}]),
        "core_edges": core_edges, "candidate_augmented_edges": edges,
        "retained_wp_groups": manifest["groups"], "retained_wp_interactions": manifest["interactions"],
        "tf_decoder": decoder, "alteration_sources": ALTERATION_SOURCES, "reachability": reachability,
        "audit": {
            "signalling_nodes": len(signalling_entities), "core_edges": len(core_edges), "candidate_augmented_edges": len(edges),
            "wp_gene_edges_after_group_expansion": sum("WikiPathways:WP4155" in x["provenance"] for x in edges),
            "augmented_edges": sum("WikiPathways:WP4155" not in x["provenance"] for x in edges),
            "core_positive_edges": sum(x["sign"] == 1 for x in core_edges), "core_negative_edges": sum(x["sign"] == -1 for x in core_edges),
            "core_neutral_edges": sum(x["sign"] == 0 for x in core_edges),
            "ambiguous_signed_pairs": sum(len(x) > 1 for x in ambiguous_pairs.values()),
            "terminal_tfs": terminal_tfs, "decoder_edges": len(decoder),
            "decoder_targets": len({x["target"] for x in decoder}),
        },
        "warnings": [
            "WikiPathways family/group edges are expanded to member genes but original group objects are retained.",
            "Neutral WP binding/conversion edges remain present but must not be assigned an activation sign.",
            "The TF decoder is not a signalling layer and must be regularised/evaluated separately.",
            "Only WP4155 edges plus three explicit ERBB2 bridge pairs form the core propagation graph; all other within-node OmniPath edges remain candidate augmentations.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"audit": payload["audit"], "reachability": reachability}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
