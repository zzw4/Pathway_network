# ARISE source-specific Gate A/B audit

> **Status correction:** This report is exploratory. Its formal pass/fail decisions are superseded by `MULTIOMICS_EXECUTION_PLAN_2026-10-02.md`: gene dependence invalidates the independent-binomial claim gate, and graph mechanism/complex semantics remain unresolved. Keep the numerical results, but do not use them as source admission decisions.

**Date:** 2 October 2026  
**Controlling documents:** `alteration_pathway_expression_project_plan.md` and `SOURCE_SPECIFIC_EXECUTION_CORRECTION_2026-10-02.md`

## Scope discipline

This audit tests the frozen source-resolved main line. It does not reconstruct global expression, pool alteration sources, use PROGENy genes as direct targets, fit a neural model, or claim causal effects.

Each source is injected separately with sign `+1`, propagated edge by edge through the signed WP4155 core graph, and decoded only through signed direct DoRothEA TF-target edges. Simple paths are enumerated up to eight edges. A TF is excluded if both positive and negative source-to-TF paths exist; a target is excluded if its retained TF evidence gives conflicting predicted signs.

## Gate A — source-specific graph manifest

The machine-readable manifest is `outputs/source_specific_graph_manifest.json`.

| Source | Functional event | Reachable signalling nodes | Reachable TFs | Directional direct targets | Sign conflicts |
|---|---|---:|---|---:|---:|
| KRAS | G12*, G13*, Q61* | 9 | ELK1, FOS | 124 | 0 |
| CTNNB1 | D32/S33/G34/S37/T41/S45 | 6 | LEF1, MYC, TCF7, TCF7L2 | 405 (403 Entrez-mapped) | 0 |

The principal explicit paths are:

- KRAS → ARAF/BRAF/RAF1 → MAP2K1/2 → MAPK1/3 → ELK1 → FOS, all positive.
- CTNNB1 → LEF1/TCF7/TCF7L2, and CTNNB1 → TCF/LEF → MYC, all positive in the frozen graph.

Locked intermediate readouts remain separate from RNA fitting:

- KRAS: ERK and MEK activation-site phosphoprotein modules.
- CTNNB1: total beta-catenin protein and mapped beta-catenin phosphosites.

**Gate A decision: pass.** The source-level graphs are explicit, directional, signed, and conflict-audited.

## Gate B — observed direction audit

TCGA-UCEC primary tumours were analysed using per-gene OLS with HC3 robust covariance. Each source coefficient was adjusted jointly for the other source, PTEN loss of function, canonical PIK3CA, PIK3R1 loss of function, damaging TP53, MSI, subtype, tumour type and grade when sufficiently observed. RNA was standardized across the analysis cohort. The sources were never pooled.

Carrier support in the expression cohort was 88 for KRAS and 100 for CTNNB1.

| Source | Tested targets | All targets matching | nominal P<0.05 matching | FDR q<0.10 matching | Median signed effect |
|---|---:|---:|---:|---:|---:|
| KRAS | 124 | 62/124 = 50.0% (P=0.536) | 15/30 = 50.0% (P=0.572) | 9/18 = 50.0% (P=0.593) | -0.0007 SD |
| CTNNB1 | 402 | 221/402 = 55.0% (P=0.0258) | 85/135 = 63.0% (P=0.00164) | 75/114 = 65.8% (P=0.000480) | +0.0260 SD |

The reported P values are one-sided exact binomial tests against 50% sign concordance. They are descriptive feasibility evidence, not a causal test and not independent external validation.

Decoder localisation gives:

- KRAS–ELK1: 17/32 overall and 4/6 at q<0.10 match.
- KRAS–FOS: 47/94 overall and 6/13 at q<0.10 match.
- CTNNB1–LEF1: 16/22 overall and 6/8 at q<0.10 match.
- CTNNB1–MYC: 197/367 overall and 64/99 at q<0.10 match.
- CTNNB1–TCF7: 5/6 overall and 4/4 at q<0.10 match.
- CTNNB1–TCF7L2: 15/21 overall and 11/13 at q<0.10 match.

## Decision

**Gate B is a source-specific split decision: CTNNB1 passes; KRAS fails.**

The positive KRAS phospho-intermediate observation cannot rescue its failed direct-target RNA direction audit, because intermediate concordance and downstream expression concordance test different required links. Conversely, the CTNNB1 result is strong enough to admit that source to Gate C, but one surviving source is insufficient for the frozen two-source project claim.

Therefore:

1. do not start the neural model;
2. admit CTNNB1 alone to the additive/linear-signed Gate C implementation;
3. keep KRAS as a documented failed source, not silently prune discordant targets;
4. before replacing KRAS, run a preregistered backup-source screen using the same untouched Gate A/B rules and existing graph; no target-level result shopping is permitted;
5. the full project proceeds only if a second source passes support, RNA direction, and non-RNA intermediate requirements.

## Reproducibility

- `scripts/build_source_specific_manifests.py`
- `scripts/audit_source_specific_directions.py`
- `outputs/source_specific_graph_manifest.json`
- `outputs/source_specific_direction_audit.json`

