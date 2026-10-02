"""Audit cBioPortal modality linkage without inspecting outcome associations."""
import argparse
import json
from pathlib import Path
import requests

API = "https://www.cbioportal.org/api"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    session = requests.Session()
    def get(path):
        response = session.get(API + path, params={"pageSize": 100000}, timeout=120)
        response.raise_for_status()
        return response.json()
    studies = [s for s in get("/studies") if "endometr" in (s.get("name", "") + s.get("description", "")).lower() or s["studyId"].startswith("ucec")]
    results = {}
    for study in studies:
        sid = study["studyId"]
        if "cptac" not in sid:
            continue
        profiles = get(f"/studies/{sid}/molecular-profiles")
        samples = get(f"/studies/{sid}/samples")
        lists = get(f"/studies/{sid}/sample-lists")
        availability = {}
        for sample_list in lists:
            lid = sample_list["sampleListId"]
            availability[lid] = get(f"/sample-lists/{lid}/sample-ids")
        mapping = {s["sampleId"]: s["patientId"] for s in samples}
        results[sid] = {
            "study": study, "sample_count": len(samples),
            "patient_count": len(set(mapping.values())),
            "profiles": profiles,
            "sample_lists": {k: {"samples": len(v), "patients": len({mapping[x] for x in v if x in mapping})} for k, v in availability.items()},
            "sample_list_ids": availability,
            "sample_to_patient": mapping,
            "clinical_attributes": get(f"/studies/{sid}/clinical-attributes"),
        }
    payload = {"candidate_studies": studies, "cptac_schema": results,
               "warning": "Sample-list membership is assay availability metadata, not proof of nonmissing feature values or cross-cohort independence."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({k: {"samples": v["sample_count"], "patients": v["patient_count"], "profiles": [p["molecularProfileId"] for p in v["profiles"]], "sample_lists": v["sample_lists"]} for k,v in results.items()}, indent=2))

if __name__ == "__main__":
    main()
