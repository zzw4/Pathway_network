#!/usr/bin/env python3
"""Convert the complete WP4155 GPML diagram into an auditable graph manifest."""

from __future__ import annotations

import argparse
import json
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path


NS = {"g": "http://pathvisio.org/GPML/2013a"}


def local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gpml", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root = ET.parse(args.gpml).getroot()

    nodes = []
    graph_ref_to_entity = {}
    group_members = defaultdict(list)
    for element in root.findall("g:DataNode", NS):
        graph_id = element.attrib["GraphId"]
        xref = element.find("g:Xref", NS)
        node = {
            "id": f"node:{graph_id}", "graph_id": graph_id,
            "label": element.attrib.get("TextLabel", "").strip(), "type": element.attrib.get("Type"),
            "group_ref": element.attrib.get("GroupRef"),
            "xref_database": xref.attrib.get("Database") if xref is not None else None,
            "xref_id": xref.attrib.get("ID") if xref is not None else None,
        }
        nodes.append(node)
        graph_ref_to_entity[graph_id] = node["id"]
        if node["group_ref"]:
            group_members[node["group_ref"]].append(node["id"])

    groups = []
    group_id_to_entity = {}
    for element in root.findall("g:Group", NS):
        group_id = element.attrib["GroupId"]
        entity_id = f"group:{group_id}"
        group = {
            "id": entity_id, "group_id": group_id, "graph_id": element.attrib.get("GraphId"),
            "style": element.attrib.get("Style", "Group"), "members": group_members.get(group_id, []),
        }
        groups.append(group)
        group_id_to_entity[group_id] = entity_id
        if group["graph_id"]:
            graph_ref_to_entity[group["graph_id"]] = entity_id

    # Interactions can themselves be branch points referenced by another interaction.
    interaction_elements = root.findall("g:Interaction", NS)
    for element in interaction_elements:
        graph_ref_to_entity[element.attrib["GraphId"]] = f"junction:{element.attrib['GraphId']}"

    interactions = []
    unresolved = Counter()
    arrowheads = Counter()
    for element in interaction_elements:
        graph_id = element.attrib["GraphId"]
        graphics = element.find("g:Graphics", NS)
        points = graphics.findall("g:Point", NS) if graphics is not None else []
        if len(points) < 2:
            unresolved["fewer_than_two_points"] += 1
            continue
        first, last = points[0], points[-1]
        first_arrow, last_arrow = first.attrib.get("ArrowHead"), last.attrib.get("ArrowHead")
        arrow = last_arrow or first_arrow
        arrowheads[arrow or "none"] += 1
        first_ref, last_ref = first.attrib.get("GraphRef"), last.attrib.get("GraphRef")
        source_ref, target_ref = (last_ref, first_ref) if first_arrow and not last_arrow else (first_ref, last_ref)
        source = graph_ref_to_entity.get(source_ref) if source_ref else None
        target = graph_ref_to_entity.get(target_ref) if target_ref else None
        if source is None:
            unresolved["source_missing_or_coordinate_only"] += 1
        if target is None:
            unresolved["target_missing_or_coordinate_only"] += 1
        directed = arrow not in {None, "mim-binding"}
        sign = -1 if arrow == "mim-inhibition" else (1 if arrow in {"Arrow", "mim-stimulation", "mim-transcription-translation"} else 0)
        interactions.append({
            "id": f"interaction:{graph_id}", "graph_id": graph_id,
            "source": source, "target": target, "source_graph_ref": source_ref, "target_graph_ref": target_ref,
            "arrowhead": arrow, "directed": directed, "sign": sign,
            "line_style": graphics.attrib.get("LineStyle") if graphics is not None else None,
            "point_count": len(points), "fully_resolved": source is not None and target is not None,
        })

    pathways = [x for x in nodes if x["type"] == "Pathway"]
    genes = [x for x in nodes if x["type"] == "GeneProduct"]
    gene_labels = sorted({x["label"] for x in genes})
    payload = {
        "pathway": {
            "id": "WP4155", "name": root.attrib.get("Name"), "version": root.attrib.get("Version"),
            "last_modified": root.attrib.get("Last-Modified"), "organism": root.attrib.get("Organism"),
        },
        "nodes": nodes, "groups": groups, "interactions": interactions,
        "membership_edges": [{"source": member, "target": group["id"], "relation": "member_of"} for group in groups for member in group["members"]],
        "audit": {
            "data_nodes": len(nodes), "gene_node_instances": len(genes), "unique_gene_symbols": len(gene_labels),
            "gene_symbols": gene_labels, "pathway_module_labels": [x["label"] for x in pathways],
            "groups": len(groups), "interactions": len(interactions),
            "fully_resolved_interactions": sum(x["fully_resolved"] for x in interactions),
            "directed_interactions": sum(x["directed"] for x in interactions),
            "signed_positive_interactions": sum(x["sign"] == 1 for x in interactions),
            "signed_negative_interactions": sum(x["sign"] == -1 for x in interactions),
            "neutral_or_unspecified_interactions": sum(x["sign"] == 0 for x in interactions),
            "arrowheads": dict(arrowheads), "unresolved_endpoint_reasons": dict(unresolved),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(payload["audit"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
