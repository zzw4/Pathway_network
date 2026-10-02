"""Fetch a scoped CPTAC development package; preserve missingness and provenance."""
import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
import requests
from audit_ucec_expression_signal import hgnc_map

API = "https://www.cbioportal.org/api"
STUDY = "ucec_cptac_2020"
# Functional annotation review is separate from availability; unknown sites are not activity labels.
SITES = ["MAPK1_T185", "MAPK1_Y187", "MAPK3_T202", "MAPK3_Y204",
         "MAP2K1_S218", "MAP2K1_S222", "MAP2K2_S222", "MAP2K2_S226",
         "AKT1_T308", "AKT1_S473", "AKT2_T309", "AKT2_S474",
         "FOXO3_T32", "FOXO3_S253", "FOXO3_S315", "GSK3B_S9",
         "CTNNB1_S33", "CTNNB1_S37", "CTNNB1_T41", "CTNNB1_S45"]

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--graph", type=Path, required=True)
    p.add_argument("--decoder", type=Path, required=True)
    p.add_argument("--cache-dir", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    a = p.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = a.out_dir / "raw_api"
    raw_dir.mkdir(exist_ok=True)
    def fetch(key, path, body=None, params=None):
        cache = raw_dir / (key + ".json")
        if cache.exists():
            return json.loads(cache.read_text())
        response = requests.get(API+path, params=params, timeout=240) if body is None else requests.post(API+path, params=params, json=body, timeout=240)
        if response.status_code == 404:
            result = {"unavailable": True, "status": 404}
        else:
            response.raise_for_status()
            result = response.json()
        cache.write_text(json.dumps(result), encoding="utf-8")
        return result
    graph = json.loads(a.graph.read_text())
    decoder = json.loads(a.decoder.read_text())
    mapping = hgnc_map(a.cache_dir / "hgnc_complete_set.txt")
    core = sorted({n["id"] for n in graph["nodes"] if n["id"] in mapping})
    targets = sorted({e["target"] for e in decoder["direct_tf_edges"] if e["sign"] in [-1,1] and e["target"] in mapping})
    samples = fetch("samples", f"/studies/{STUDY}/samples", params={"pageSize":100000})
    sample_ids = sorted(x["sampleId"] for x in samples)
    patient_map = {x["sampleId"]:x["patientId"] for x in samples}
    if len(set(patient_map.values())) != len(sample_ids):
        raise ValueError("Multiple samples per patient: explicit selection required")
    sequenced = set(fetch("sequenced", f"/sample-lists/{STUDY}_sequenced/sample-ids"))
    mutations = fetch("mutations", f"/molecular-profiles/{STUDY}_mutations/mutations/fetch", {"sampleIds":sample_ids,"entrezGeneIds":[mapping[x] for x in core]}, {"projection":"DETAILED"})
    attrs = fetch("clinical_attributes", f"/studies/{STUDY}/clinical-attributes")
    sample_clin = fetch("sample_clinical", f"/studies/{STUDY}/clinical-data", params={"clinicalDataType":"SAMPLE","pageSize":100000})
    patient_clin = fetch("patient_clinical", f"/studies/{STUDY}/clinical-data", params={"clinicalDataType":"PATIENT","pageSize":100000})
    panels = {"rna": ("mrna",sorted(set(core+targets))), "protein": ("protein_quantification",core), "cnv": ("gistic",core)}
    arrays, report = {}, {}
    for modality,(suffix, genes) in panels.items():
        records=[]
        for start in range(0,len(genes),150):
            records.extend(fetch(f"{modality}_{start}",f"/molecular-profiles/{STUDY}_{suffix}/molecular-data/fetch", {"sampleIds":sample_ids,"entrezGeneIds":[mapping[x] for x in genes[start:start+150]]}))
        frame=pd.DataFrame(records)
        if frame.empty:
            matrix=pd.DataFrame(index=sample_ids,columns=[mapping[x] for x in genes],dtype=float)
        else:
            frame["value"]=pd.to_numeric(frame["value"],errors="coerce")
            matrix=frame.pivot_table(index="sampleId",columns="entrezGeneId",values="value",aggfunc="first").reindex(index=sample_ids,columns=[mapping[x] for x in genes])
        arrays[modality]=matrix.to_numpy(dtype=float)
        arrays[modality+"_genes"]=np.array(genes)
        report[modality]={"shape":list(matrix.shape),"patients_with_any":int(matrix.notna().any(axis=1).sum()),"observed_fraction":float(matrix.notna().to_numpy().mean()),"feature_nonmissing":dict(zip(genes,matrix.notna().sum().astype(int).tolist()))}
    def site_fetch(site):
        return site,fetch("phospho_"+site,f"/generic-assay-data/{STUDY}_phosphoproteome/generic-assay/{site}")
    phospho=pd.DataFrame(index=sample_ids,columns=SITES,dtype=float)
    absent=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for site,data in pool.map(site_fetch,SITES):
            if isinstance(data,dict):
                absent.append(site)
                continue
            values={r["sampleId"]:pd.to_numeric(r.get("value"),errors="coerce") for r in data}
            phospho[site]=pd.Series(values).reindex(sample_ids)
    arrays["phospho"]=phospho.to_numpy(dtype=float)
    arrays["phospho_sites"]=np.array(SITES)
    arrays["sample_ids"]=np.array(sample_ids)
    arrays["patient_ids"]=np.array([patient_map[s] for s in sample_ids])
    arrays["sequenced"]=np.array([s in sequenced for s in sample_ids])
    report["phospho"]={"shape":list(phospho.shape),"patients_with_any":int(phospho.notna().any(axis=1).sum()),"feature_nonmissing":phospho.notna().sum().astype(int).to_dict(),"absent_sites":absent}
    joint=np.logical_and.reduce([np.isfinite(arrays[k]).any(axis=1) for k in panels]+[phospho.notna().any(axis=1).to_numpy(),arrays["sequenced"]])
    dataset=a.out_dir/"paired_development.npz"
    np.savez_compressed(dataset,**arrays)
    report.update({"study":STUDY,"patients":len(sample_ids),"patients_with_any_in_all_modalities":int(joint.sum()),"core_genes":len(core),"direct_target_genes":len(targets),"mutation_records":len(mutations),"sha256":hashlib.sha256(dataset.read_bytes()).hexdigest(),"clinical_attribute_definitions":attrs,"policy":"Raw values only; no global normalisation, imputation, association testing or train/test splitting."})
    (a.out_dir/"coverage_report.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k not in ["clinical_attribute_definitions","rna","protein","cnv"]},indent=2))

if __name__ == "__main__":
    main()
