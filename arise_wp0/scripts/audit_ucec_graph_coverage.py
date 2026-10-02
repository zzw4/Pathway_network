#!/usr/bin/env python3
"""Audit whether a signed UCEC signalling graph connects alteration sources to RNA targets."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict, deque
from pathlib import Path

import pandas as pd


SOURCES = ["KRAS", "CTNNB1", "TP53", "PIK3R1", "PTEN", "PIK3CA"]
PATHWAYS = ["MAPK", "WNT", "p53", "PI3K"]


def clean_symbol(value: object) -> str | None:
    if not isinstance(value, str) or not value or "_" in value or ":" in value:
        return None
    return value


def bfs(adjacency: dict[str, set[str]], source: str, max_depth: int) -> dict[str, int]:
    distance = {source: 0}
    queue = deque([source])
    while queue:
        node = queue.popleft()
        if distance[node] >= max_depth:
            continue
        for nxt in adjacency.get(node, set()):
            if nxt not in distance:
                distance[nxt] = distance[node] + 1
                queue.append(nxt)
    return distance


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--omnipath", required=True, type=Path)
    parser.add_argument("--dorothea", required=True, type=Path)
    parser.add_argument("--progeny", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--max-depth", type=int, default=6)
    args = parser.parse_args()

    signalling = pd.read_csv(args.omnipath, sep="\t", low_memory=False)
    signalling = signalling[(signalling["is_directed"] == True) & ((signalling["is_stimulation"] == True) | (signalling["is_inhibition"] == True))].copy()  # noqa: E712
    signalling["src"] = signalling["source_genesymbol"].map(clean_symbol)
    signalling["dst"] = signalling["target_genesymbol"].map(clean_symbol)
    signalling = signalling.dropna(subset=["src", "dst"])
    adjacency = defaultdict(set)
    edge_signs = defaultdict(set)
    for row in signalling.itertuples():
        adjacency[row.src].add(row.dst)
        if row.is_stimulation:
            edge_signs[(row.src, row.dst)].add(1)
        if row.is_inhibition:
            edge_signs[(row.src, row.dst)].add(-1)

    regulons = pd.read_csv(args.dorothea, sep="\t", low_memory=False)
    regulons = regulons[regulons["dorothea_level"].astype(str).str.contains("A|B", regex=True)].copy()
    regulons["tf"] = regulons["source_genesymbol"].map(clean_symbol)
    regulons["target_gene"] = regulons["target_genesymbol"].map(clean_symbol)
    regulons = regulons.dropna(subset=["tf", "target_gene"])
    tf_set = set(regulons["tf"])

    raw = pd.read_csv(args.progeny, sep="\t")
    model = raw.pivot_table(index=["genesymbol", "record_id"], columns="label", values="value", aggfunc="first").reset_index()
    model["p_value"] = pd.to_numeric(model["p_value"], errors="coerce")
    targets = {}
    for pathway in PATHWAYS:
        targets[pathway] = set(model[model["pathway"] == pathway].sort_values("p_value").drop_duplicates("genesymbol").head(100)["genesymbol"])

    source_results = {}
    union_reachable_tfs = set()
    for source in SOURCES:
        distances = bfs(adjacency, source, args.max_depth)
        reachable_tfs = tf_set & set(distances)
        union_reachable_tfs |= reachable_tfs
        regulated = set(regulons.loc[regulons["tf"].isin(reachable_tfs), "target_gene"])
        source_results[source] = {
            "reachable_signalling_nodes": len(distances), "reachable_tfs": len(reachable_tfs),
            "nearest_tfs": sorted([{"tf": tf, "distance": distances[tf]} for tf in reachable_tfs], key=lambda x: (x["distance"], x["tf"]))[:30],
            "progeny_target_coverage": {p: {"covered": len(t & regulated), "total": len(t), "fraction": len(t & regulated) / len(t)} for p, t in targets.items()},
        }

    union_regulated = set(regulons.loc[regulons["tf"].isin(union_reachable_tfs), "target_gene"])
    payload = {
        "sources": SOURCES, "max_signalling_depth": args.max_depth,
        "graph": {"signed_directed_edges": int(len(signalling)), "nodes": len(set(signalling["src"]) | set(signalling["dst"])),
                  "ambiguous_sign_edges": sum(len(v) > 1 for v in edge_signs.values()),
                  "dorothea_ab_edges": int(len(regulons)), "dorothea_tfs": len(tf_set)},
        "source_results": source_results,
        "union": {"reachable_tfs": len(union_reachable_tfs),
                  "progeny_target_coverage": {p: {"covered": len(t & union_regulated), "total": len(t), "fraction": len(t & union_regulated) / len(t)} for p, t in targets.items()}},
        "gate": "Each validated source must reach TFs, and the union must cover >=50% of top-100 targets in at least three pathway programmes.",
    }
    validated = ["KRAS", "CTNNB1"]
    payload["passes_gate"] = all(source_results[x]["reachable_tfs"] > 0 for x in validated) and sum(v["fraction"] >= 0.5 for v in payload["union"]["progeny_target_coverage"].values()) >= 3
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
