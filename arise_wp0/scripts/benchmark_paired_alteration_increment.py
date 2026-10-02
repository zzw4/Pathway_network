"""Patient-level nested ridge pilot; establishes signal, not pathway causality."""
import argparse
import json
import re
from pathlib import Path
import numpy as np
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from refine_functional_sources import matches, is_lof, TP53_HOTSPOTS

FEATURES=["KRAS_CANONICAL","CTNNB1_EXON3","PIK3CA_CANONICAL","PTEN_LOF","PIK3R1_LOF","TP53_LOF","TP53_HOTSPOT"]

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--data-dir",type=Path,required=True)
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    d=np.load(a.data_dir/"paired_development.npz",allow_pickle=False)
    patients=d["patient_ids"].astype(str)
    samples=d["sample_ids"].astype(str)
    index={s:i for i,s in enumerate(samples)}
    X=np.zeros((len(samples),len(FEATURES)))
    records=json.loads((a.data_dir/"raw_api/mutations.json").read_text())
    for r in records:
        i=index.get(r.get("sampleId"))
        if i is None:continue
        gene=(r.get("gene") or {}).get("hugoGeneSymbol","")
        change=(r.get("proteinChange") or "").replace("p.","").upper()
        residue=re.fullmatch(r"[A-Z](\d+)[A-Z]",change)
        for j,f in enumerate(FEATURES):
            hit=matches(f,gene,r)
            if f=="CTNNB1_EXON3":hit=gene=="CTNNB1" and bool(residue and int(residue.group(1)) in [32,33,34,37,41,45])
            if f=="TP53_LOF":hit=gene=="TP53" and is_lof(r)
            if f=="TP53_HOTSPOT":hit=gene=="TP53" and change in TP53_HOTSPOTS
            if hit:X[i,j]=1
    if not d["sequenced"].all():raise ValueError("Unsequenced patients require exclusion")
    # Purity is a sensitivity covariate, not a claim of causal adjustment.
    clinical=json.loads((a.data_dir/"raw_api/sample_clinical.json").read_text())
    values={}
    for r in clinical:
        if r["clinicalAttributeId"]=="PURITY_CANCER":
            try:values[r["sampleId"]]=float(r["value"])
            except (ValueError,TypeError):pass
    C=np.array([[values.get(s,np.nan)] for s in samples])
    # No output-derived subtype, proliferation signature or observed mediator input.
    manifest=json.loads(a.manifest.read_text())
    rna_genes=d["rna_genes"].astype(str).tolist()
    panels={s:[rna_genes.index(t["target"]) for t in manifest["sources"][s]["directional_targets"] if t["target"] in rna_genes and t["target"]!=s] for s in ["KRAS","CTNNB1"]}
    panels.update({"protein_CTNNB1":[d["protein_genes"].astype(str).tolist().index("CTNNB1")],"phospho_ERK":[d["phospho_sites"].astype(str).tolist().index(s) for s in ["MAPK1_T185","MAPK1_Y187","MAPK3_T202","MAPK3_Y204"]],"phospho_FOXO3":[d["phospho_sites"].astype(str).tolist().index(s) for s in ["FOXO3_T32","FOXO3_S253"]]})
    alphas=[1.,10.,100.]
    outer=list(KFold(5,shuffle=True,random_state=20261002).split(X))
    def prep(xx,train,test):
        med=np.nanmedian(xx[train],axis=0)
        med=np.where(np.isfinite(med),med,0.)
        a1=np.where(np.isfinite(xx[train]),xx[train],med)
        a2=np.where(np.isfinite(xx[test]),xx[test],med)
        scale=StandardScaler().fit(a1)
        return scale.transform(a1),scale.transform(a2)
    results={}
    for panel,cols in panels.items():
        modality="rna" if panel in ["KRAS","CTNNB1"] else "protein" if panel.startswith("protein") else "phospho"
        Y=d[modality][:,cols]
        # Complete outcomes only for a simple auditable pilot; separate missing-site models later.
        keep=np.isfinite(Y).all(axis=0)
        Y=Y[:,keep]
        if not Y.shape[1]:
            results[panel]={"status":"no outcome complete in all patients; excluded from complete-case pilot"};continue
        model_inputs=[("purity",C),("purity_alterations",np.column_stack([C,X]))]
        source_feature={"KRAS":"KRAS_CANONICAL","CTNNB1":"CTNNB1_EXON3","protein_CTNNB1":"CTNNB1_EXON3"}.get(panel)
        if source_feature:
            remove_index=FEATURES.index(source_feature)
            model_inputs.append(("purity_other_alterations",np.column_stack([C,np.delete(X,remove_index,axis=1)])))
        predictions={name:np.zeros_like(Y) for name in ["mean"]+[name for name,_ in model_inputs]}
        null=np.zeros_like(Y)
        folds=[]
        for fold,(train,test) in enumerate(outer):
            center=Y[train].mean(axis=0);std=Y[train].std(axis=0)
            std=np.where(std>1e-8,std,1.)
            yt=(Y[train]-center)/std
            null[test]=center
            predictions["mean"][test]=center
            for name,xx in model_inputs:
                scores=[]
                for alpha in alphas:
                    errors=[]
                    for it,iv in KFold(3,shuffle=True,random_state=fold).split(train):
                        tr,va=train[it],train[iv]
                        xt,xv=prep(xx,tr,va)
                        mu=Y[tr].mean(axis=0);sd=Y[tr].std(axis=0);sd=np.where(sd>1e-8,sd,1.)
                        fit=Ridge(alpha=alpha).fit(xt,(Y[tr]-mu)/sd)
                        predicted=np.asarray(fit.predict(xv)).reshape(len(va),Y.shape[1])
                        errors.append(np.mean(((predicted-(Y[va]-mu)/sd))**2))
                    scores.append(np.mean(errors))
                alpha=alphas[int(np.argmin(scores))]
                xt,xv=prep(xx,train,test)
                fit=Ridge(alpha=alpha).fit(xt,yt)
                predicted=np.asarray(fit.predict(xv)).reshape(len(test),Y.shape[1])
                predictions[name][test]=predicted*std+center
                folds.append({"fold":fold,"model":name,"alpha":alpha})
        # Normalise squared errors using per-fold training outcome variance.
        scale=np.zeros_like(Y)
        for tr,te in outer:
            sd=Y[tr].std(axis=0);scale[te]=np.where(sd>1e-8,sd,1.)
        den=((Y-null)/scale)**2
        errors={k:((Y-v)/scale)**2 for k,v in predictions.items()}
        r2={k:float(1-v.mean()/den.mean()) for k,v in errors.items()}
        delta_patient=(errors["purity"]-errors["purity_alterations"]).mean(axis=1)
        rng=np.random.default_rng(20261002)
        boot=np.array([delta_patient[rng.integers(0,len(X),len(X))].mean()/den.mean() for _ in range(2000)])
        results[panel]={"complete_features":int(Y.shape[1]),"requested_features":len(cols),"patients":len(X),"oof_r2":r2,"delta_r2_vs_purity":r2["purity_alterations"]-r2["purity"],"patient_bootstrap_ci95":np.quantile(boot,[.025,.975]).tolist(),"selected_alphas":folds}
        if source_feature:
            carriers=X[:,FEATURES.index(source_feature)]==1
            improvement=(errors["purity_other_alterations"]-errors["purity_alterations"]).mean(axis=1)
            boot=np.array([improvement[rng.integers(0,len(X),len(X))].mean()/den.mean() for _ in range(2000)])
            results[panel]["source_feature"]=source_feature
            results[panel]["source_increment_r2"]=r2["purity_alterations"]-r2["purity_other_alterations"]
            results[panel]["source_increment_bootstrap_ci95"]=np.quantile(boot,[.025,.975]).tolist()
            results[panel]["carrier_mean_standardised_mse_improvement"]=float(improvement[carriers].mean())
    counts={f:int(X[:,j].sum()) for j,f in enumerate(FEATURES)}
    np.savez_compressed(a.data_dir/"functional_source_features.npz",X=X,feature_names=np.array(FEATURES),patient_ids=patients,sample_ids=samples,purity=C)
    payload={"carrier_counts":counts,"results":results,"policy":"Development nested patient CV pilot. Bootstrap conditions on fitted folds. Purity-only control is incomplete confounding control. Not a causal or graph/neural performance test. No CNV attribution until focal-region audit."}
    a.output.write_text(json.dumps(payload,indent=2))
    print(json.dumps({"carrier_counts":counts,"results":{k:{x:y for x,y in v.items() if x!="selected_alphas"} for k,v in results.items()}},indent=2))

if __name__=="__main__":main()
