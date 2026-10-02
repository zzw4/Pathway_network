"""Flag graph relations requiring state or complex semantics before modelling."""
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--graph",type=Path,required=True)
    p.add_argument("--wp",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    graph=json.loads(a.graph.read_text())
    wp=json.loads(a.wp.read_text())
    nodes={n["id"]:n for n in wp["nodes"]}
    groups={g["id"]:g for g in wp["groups"]}
    def resolve(entity):
        if entity in nodes:
            return [nodes[entity]["label"]]
        return [nodes[x]["label"] for x in groups.get(entity,{}).get("members",[]) if x in nodes]
    support=defaultdict(list)
    for interaction in wp["interactions"]:
        complex_endpoint=any(groups.get(interaction.get(k),{}).get("style")=="Complex" for k in ["source","target"])
        for u in resolve(interaction.get("source")):
            for v in resolve(interaction.get("target")):
                support[(u,v)].append({"interaction_id":interaction["id"],"complex_endpoint":complex_endpoint,"arrowhead":interaction["arrowhead"]})
    signs=defaultdict(set)
    for e in graph["candidate_augmented_edges"]:
        if e["sign"] in [-1,1]:
            signs[(e["source"],e["target"])].add(e["sign"])
    rows=[]
    for e in graph["core_edges"]:
        pair=(e["source"],e["target"])
        evid=support[pair]
        flags=[]
        if not e["directed"] or e["sign"]==0:
            flags.append("neutral_or_undirected")
        if len(signs[pair])>1:
            flags.append("sign_conflict_in_candidate_evidence")
        if any(x["complex_endpoint"] for x in evid):
            flags.append("complex_expanded_to_member_edges")
        arrow=e.get("wp_arrowhead") or ""
        if "transcription" in arrow:
            state="regulator_activity_to_target_abundance"
        elif e["target"]=="PIP3":
            state="enzyme_activity_to_metabolite_state"
        elif e["source"]=="PIP3":
            state="metabolite_state_to_protein_activity"
        else:
            state="unresolved_activity_abundance_or_complex_relation"
            flags.append("mechanism_and_states_require_curation")
        if e["source"]=="CTNNB1" and e["target"] in ["TCF7","TCF7L1","TCF7L2","LEF1"]:
            flags.append("beta_catenin_TCF_complex_requires_explicit_representation")
        rows.append({**e,"proposed_state_relation":state,"flags":flags,"wp_interaction_support":evid,"ready_for_activity_training":False})
    result={"policy":"Automated flags are curation candidates, not validated biochemical assignments.","core_edges":rows,"summary":{"core_edges":len(rows),"flags":dict(Counter(f for r in rows for f in r["flags"]))}}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result["summary"],indent=2))

if __name__=="__main__":
    main()
