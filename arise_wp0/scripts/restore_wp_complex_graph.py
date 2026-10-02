"""Restore original WP entities; do not treat complex membership as activation."""
import argparse
import json
from pathlib import Path
from collections import Counter, defaultdict, deque

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--wp",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    wp=json.loads(a.wp.read_text())
    entities={}
    for n in wp["nodes"]:
        if n["type"]=="GeneProduct" or n["label"]=="PIP3":
            entities[n["id"]]={"id":n["id"],"label":n["label"],"type":n["type"],"state_semantics":"requires_curation"}
    for g in wp["groups"]:
        members=[m for m in g["members"] if m in entities]
        if members:
            entities[g["id"]]={"id":g["id"],"type":"complex" if g["style"]=="Complex" else "family_or_group", "members":members,"label":" / ".join(entities[m]["label"] for m in members),"state_semantics":"requires_curation"}
    relations=[]
    membership=[]
    for e in entities.values():
        for m in e.get("members",[]):
            membership.append({"member":m,"container":e["id"],"relation":"membership","sign":None,"propagation":False})
    for e in wp["interactions"]:
        if e.get("source") not in entities or e.get("target") not in entities:
            continue
        arrow=e.get("arrowhead","") or ""
        category="transcription_to_abundance" if "transcription" in arrow else "binding_or_unsigned" if not e["directed"] or not e["sign"] else "signed_relation_mechanism_unresolved"
        relations.append({**e,"category":category,"source_label":entities[e["source"]]["label"],"target_label":entities[e["target"]]["label"],"ready_for_training":False})
    # Structural walk includes entry to a family/complex ONLY as uncertain membership.
    adj=defaultdict(list)
    for e in relations:
        if e["directed"] and e["sign"] in [-1,1]:
            adj[e["source"]].append((e["target"],False))
    for e in membership:
        adj[e["member"]].append((e["container"],True))
        adj[e["container"]].append((e["member"],True))
    reach={}
    for source in ["KRAS","CTNNB1","PTEN","PIK3CA","PIK3R1","TP53"]:
        starts=[n["id"] for n in entities.values() if n.get("label")==source]
        visited=set(starts)
        q=deque(starts)
        while q:
            u=q.popleft()
            for v,_ in adj[u]:
                if v not in visited:
                    visited.add(v);q.append(v)
        reach[source]={"structural_reachable_entities":len(visited),"labels":sorted(set(entities[x]["label"] for x in visited)),"warning":"Includes unsigned membership; this is NOT signed biological propagation reachability."}
    result={"entities":list(entities.values()),"relations":relations,"membership":membership,"structural_reachability":reach,"summary":{"entities":len(entities),"complexes":sum(e["type"]=="complex" for e in entities.values()),"families_or_groups":sum(e["type"]=="family_or_group" for e in entities.values()),"relations":len(relations),"categories":dict(Counter(e["category"] for e in relations))},"policy":"Original group entities retained. Membership is neither an activation edge nor independent member-level causal evidence. Structural reachability requires biochemical curation before training."}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result["summary"],indent=2))

if __name__=="__main__":main()
