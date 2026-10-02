#!/usr/bin/env python3
"""Train the minimal signed sparse WP4155 propagation prototype and random-graph control."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import torch
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold


MODULE_TERMINALS = {
    "MAPK": {"ELK1": 0.5, "FOS": 0.5},
    "WNT": {"LEF1": 0.25, "MYC": 0.25, "TCF7": 0.25, "TCF7L2": 0.25},
    "PI3K": {"FOXO3": -1.0},
    "p53": {"TP53": 1.0},
    "EGFR": {"EGFR": 0.5, "ERBB2": 0.5},
}
MUTATION_INJECTION = {
    "KRAS_CANONICAL": ("KRAS", 1.0), "CTNNB1_ACTIVATING": ("CTNNB1", 1.0),
    "TP53_DAMAGING": ("TP53", -1.0), "PIK3R1_LOF": ("PIK3R1", -1.0),
    "PTEN_LOF": ("PTEN", -1.0), "PIK3CA_CANONICAL": ("PIK3CA", 1.0),
}


class SignedPropagation(torch.nn.Module):
    def __init__(self, clinical_dim: int, source_injection: torch.Tensor, edge_src: torch.Tensor, edge_dst: torch.Tensor,
                 edge_sign: torch.Tensor, module_matrix: torch.Tensor, decoder_src: torch.Tensor,
                 decoder_dst: torch.Tensor, decoder_sign: torch.Tensor, n_outputs: int, steps: int):
        super().__init__()
        self.register_buffer("source_injection", source_injection)
        self.register_buffer("edge_src", edge_src)
        self.register_buffer("edge_dst", edge_dst)
        self.register_buffer("edge_sign", edge_sign)
        self.register_buffer("module_matrix", module_matrix)
        self.register_buffer("decoder_src", decoder_src)
        self.register_buffer("decoder_dst", decoder_dst)
        self.register_buffer("decoder_sign", decoder_sign)
        self.edge_raw = torch.nn.Parameter(torch.full((len(edge_src),), -1.0))
        self.decoder_raw = torch.nn.Parameter(torch.full((len(decoder_src),), -1.5))
        self.source_scale_raw = torch.nn.Parameter(torch.zeros(source_injection.shape[0]))
        self.bio_gate_raw = torch.nn.Parameter(torch.tensor(-4.0))
        self.self_raw = torch.nn.Parameter(torch.tensor(0.0))
        self.clinical = torch.nn.Linear(clinical_dim, n_outputs)
        torch.nn.init.zeros_(self.clinical.weight)
        torch.nn.init.zeros_(self.clinical.bias)
        self.bias = torch.nn.Parameter(torch.zeros(n_outputs))
        self.steps = steps
        self.n_nodes = source_injection.shape[1]

    def forward(self, clinical: torch.Tensor, sources: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        scale = torch.nn.functional.softplus(self.source_scale_raw)
        h0 = (sources * scale) @ self.source_injection
        h = h0
        weights = self.edge_sign * torch.nn.functional.softplus(self.edge_raw)
        self_weight = torch.sigmoid(self.self_raw)
        for _ in range(self.steps):
            messages = h[:, self.edge_src] * weights
            aggregate = torch.zeros((h.shape[0], self.n_nodes), device=h.device, dtype=h.dtype)
            aggregate.index_add_(1, self.edge_dst, messages)
            h = torch.tanh(self_weight * h + aggregate + h0)
        modules = h @ self.module_matrix.T
        decoder_sources = torch.cat([h, modules], dim=1)
        dweights = self.decoder_sign * torch.nn.functional.softplus(self.decoder_raw)
        contributions = decoder_sources[:, self.decoder_src] * dweights
        biological = torch.zeros((h.shape[0], self.bias.numel()), device=h.device, dtype=h.dtype)
        biological.index_add_(1, self.decoder_dst, contributions)
        prediction = self.clinical(clinical) + torch.sigmoid(self.bio_gate_raw) * biological + self.bias
        penalty = torch.mean(torch.abs(weights)) + torch.mean(torch.abs(dweights))
        return prediction, penalty


def build_tensors(graph: dict, layer: dict, genes: list[str], source_names: list[str], randomise: bool, seed: int):
    nodes = [x["id"] for x in graph["nodes"]]
    node_index = {x: i for i, x in enumerate(nodes)}
    output_index = {x: i for i, x in enumerate(genes)}
    source_injection = np.zeros((len(source_names), len(nodes)), dtype=np.float32)
    for i, name in enumerate(source_names):
        if name in MUTATION_INJECTION:
            gene, value = MUTATION_INJECTION[name]
        elif name.startswith("CNA_"):
            gene, value = name[4:], 1.0
        else:
            continue
        if gene in node_index:
            source_injection[i, node_index[gene]] = value

    edge_rows = [x for x in graph["core_edges"] if x["directed"] and x["sign"] != 0 and x["source"] in node_index and x["target"] in node_index]
    edge_src = np.array([node_index[x["source"]] for x in edge_rows], dtype=np.int64)
    edge_dst = np.array([node_index[x["target"]] for x in edge_rows], dtype=np.int64)
    if randomise:
        permutation = np.random.default_rng(seed).permutation(len(nodes))
        edge_src, edge_dst = permutation[edge_src], permutation[edge_dst]
    edge_sign = np.array([x["sign"] for x in edge_rows], dtype=np.float32)

    modules = list(MODULE_TERMINALS)
    module_matrix = np.zeros((len(modules), len(nodes)), dtype=np.float32)
    for i, module in enumerate(modules):
        for gene, weight in MODULE_TERMINALS[module].items():
            if gene in node_index:
                module_matrix[i, node_index[gene]] = weight

    decoder = []
    # Direct TF edges only for the fixed primary outputs; omit unsigned direct edges.
    for x in layer["direct_tf_edges"]:
        if x["source"] in node_index and x["target"] in output_index and x["sign"] != 0:
            decoder.append((node_index[x["source"]], output_index[x["target"]], x["sign"]))
    for x in layer["footprint_edges"]:
        module = x["source"].split(":", 1)[1]
        if module in modules and x["target"] in output_index and x["sign"] != 0:
            decoder.append((len(nodes) + modules.index(module), output_index[x["target"]], x["sign"]))
    # Deduplicate exact typed/sign edges.
    decoder = sorted(set(decoder))
    return (
        torch.tensor(source_injection), torch.tensor(edge_src), torch.tensor(edge_dst), torch.tensor(edge_sign),
        torch.tensor(module_matrix), torch.tensor([x[0] for x in decoder], dtype=torch.long),
        torch.tensor([x[1] for x in decoder], dtype=torch.long), torch.tensor([x[2] for x in decoder], dtype=torch.float32),
        {"nodes": len(nodes), "core_signed_edges": len(edge_rows), "decoder_edges": len(decoder), "modules": modules},
    )


def run_condition(data, graph, layer, condition: str, args, device) -> dict:
    clinical = data["clinical"].astype(np.float32)
    sources = data["sources"].astype(np.float32)
    y = data["y"].astype(np.float32)
    base_prediction = data["clinical_oof"].astype(np.float32) if args.residualize_clinical else np.zeros_like(y)
    training_target = y - base_prediction
    if args.residualize_clinical:
        clinical = np.zeros((len(y), 1), dtype=np.float32)
    genes = data["genes"].astype(str).tolist()
    source_names = data["source_names"].astype(str).tolist()
    predictions = np.full_like(y, np.nan)
    fold_results = []
    folds = KFold(n_splits=5, shuffle=True, random_state=args.seed)
    tensors = build_tensors(graph, layer, genes, source_names, condition == "random_graph", args.seed)
    *model_tensors, audit = tensors
    model_tensors = [x.to(device) for x in model_tensors]
    for fold, (train_all, test) in enumerate(folds.split(y)):
        rng = np.random.default_rng(args.seed + fold)
        shuffled = rng.permutation(train_all)
        n_val = max(40, int(0.15 * len(shuffled)))
        val, train = shuffled[:n_val], shuffled[n_val:]
        cmean, cstd = clinical[train].mean(0), clinical[train].std(0)
        smean, sstd = sources[train].mean(0), sources[train].std(0)
        cstd[cstd == 0] = 1; sstd[sstd == 0] = 1
        c = torch.tensor((clinical - cmean) / cstd, device=device)
        s = torch.tensor((sources - smean) / sstd, device=device)
        target = torch.tensor(training_target, device=device)
        torch.manual_seed(args.seed + fold)
        model = SignedPropagation(clinical.shape[1], *model_tensors, y.shape[1], args.steps).to(device)
        optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
        best_state, best_val, patience = None, float("inf"), 0
        for epoch in range(args.epochs):
            model.train(); optimizer.zero_grad()
            pred, penalty = model(c[train], s[train])
            loss = torch.mean((pred - target[train]) ** 2) + args.l1 * penalty
            loss.backward(); optimizer.step()
            model.eval()
            with torch.no_grad():
                val_pred, _ = model(c[val], s[val])
                val_loss = float(torch.mean((val_pred - target[val]) ** 2))
            if val_loss < best_val - 1e-5:
                best_val, patience = val_loss, 0
                best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            else:
                patience += 1
                if patience >= args.patience:
                    break
        model.load_state_dict(best_state)
        model.eval()
        with torch.no_grad():
            test_pred, _ = model(c[test], s[test])
        predictions[test] = base_prediction[test] + test_pred.cpu().numpy()
        fold_results.append({"fold": fold, "epochs": epoch + 1, "best_val_mse": best_val})
    r2 = np.array([r2_score(y[:, i], predictions[:, i]) for i in range(y.shape[1])])
    return {"condition": condition, "mean_r2": float(r2.mean()), "median_r2": float(np.median(r2)),
            "genes_positive_r2": int((r2 > 0).sum()), "per_gene_r2": {g: float(v) for g, v in zip(genes, r2)},
            "folds": fold_results, "architecture": audit}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--graph", required=True, type=Path)
    parser.add_argument("--output-layer", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--epochs", type=int, default=600)
    parser.add_argument("--patience", type=int, default=50)
    parser.add_argument("--steps", type=int, default=8)
    parser.add_argument("--lr", type=float, default=0.01)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--l1", type=float, default=1e-4)
    parser.add_argument("--seed", type=int, default=20261001)
    parser.add_argument("--residualize-clinical", action="store_true")
    args = parser.parse_args()
    data = np.load(args.data, allow_pickle=False)
    graph = json.loads(args.graph.read_text(encoding="utf-8"))
    layer = json.loads(args.output_layer.read_text(encoding="utf-8"))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    results = [run_condition(data, graph, layer, x, args, device) for x in ["wp4155", "random_graph"]]
    payload = {"device": str(device), "config": vars(args) | {"data": str(args.data), "graph": str(args.graph), "output_layer": str(args.output_layer), "output": str(args.output)},
               "results": results, "passes_structure_test": results[0]["mean_r2"] > results[1]["mean_r2"]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")
    print(json.dumps({"device": str(device), "results": [{k: v for k, v in x.items() if k not in ["per_gene_r2", "folds"]} for x in results],
                      "passes_structure_test": payload["passes_structure_test"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
