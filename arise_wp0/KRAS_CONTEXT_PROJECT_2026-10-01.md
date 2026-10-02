# KRAS context-conditioned signalling propagation project

## Working question

Why does the same canonical KRAS alteration produce a conserved ERK activation core but different downstream transcriptional programmes across tumour contexts?

The computational task is not to classify KRAS status or reproduce a pathway score. It is to learn a constrained alteration-to-signalling-to-expression map with shared pathway parameters and cancer-context gates, then test whether the learned paths recover independently measured phosphoprotein intermediates and transfer to a held-out cancer.

## Empirical basis

### TCGA RNA footprint

After adjustment for available tumour type, molecular subtype and MSI covariates, canonical KRAS hotspots associate with increased PROGENy MAPK activity in all three development cohorts:

| Cohort | Complete cases | KRAS carriers | Adjusted beta | P value |
|---|---:|---:|---:|---:|
| UCEC | 505 | 88 | 1.30 | 8.7e-5 |
| LUAD | 508 | 163 | 1.46 | 2.1e-8 |
| COADREAD | 555 | 190 | 1.27 | 2.5e-7 |

The gene-level effects are only moderately conserved. Pairwise correlations across the same 100 MAPK footprint genes are 0.43–0.58. Only 8–12 genes are significant at FDR 5% in both members of a cancer pair, whereas 35–39 are significant in at least one. This supports a shared core plus context-specific output rather than a universal fixed signature.

### CPTAC phosphoprotein gate

A frozen four-site ERK activation-loop panel was adjusted for corresponding total protein and available cohort covariates.

| Cohort | KRAS carriers | ERK module beta | P value | Gate |
|---|---:|---:|---:|---|
| UCEC | 29 | 0.45 | 0.015 | pass |
| LUAD | 33 | 0.51 | 0.020 | pass |
| COAD | 27 | 0.05 | 0.818 | fail |

UCEC and LUAD each also contain individual ERK activation features passing within-cohort directional FDR 10%. COAD has a strong RNA footprint but no matching steady-state ERK phosphosite effect in this dataset. This discordance is a result to explain or validate, not something to hide.

## Claim boundary

The primary claim must be **context-conditioned alteration-associated propagation**, not proof that KRAS causally changes each inferred target in human tumours. Perturbation data may be used later as an out-of-domain validation set, but the tumour model will not be trained on perturbation responses.

## Proposed model

### Inputs

- canonical KRAS hotspot indicator and allele class;
- focal CNA and selected co-alterations as nuisance/context variables;
- tumour-context representation learned from expression features that exclude the evaluated MAPK output genes;
- clinical/technical covariates used for residualisation or adversarial control, not as predictive shortcuts.

### Architecture

1. A shared signed, directed RAS–RAF–MEK–ERK propagation backbone.
2. Low-rank context gates on edges or node transitions; shared weights remain the default and context-specific deviations are penalised.
3. A transcriptional decoder from terminal signalling/TF nodes to gene-expression residuals.
4. Optional phosphosite observation heads used only in CPTAC fine evaluation or auxiliary training experiments that are explicitly separated from the main TCGA comparison.

The smallest viable implementation should be used first: two to four propagation layers, low-rank gates, and a linear sparse decoder. A large generic GNN is not justified.

## Evaluation

### Primary tests

1. **Leave-one-cancer-out transfer:** train shared propagation on two cancers and predict KRAS-associated expression effects in the third.
2. **Within-cancer held-out prediction:** nested cross-validation with all samples from a patient confined to one fold.
3. **Intermediate recovery:** without selecting paths from phosphoproteomics, test whether node activations align with CPTAC ERK phosphosite modules.
4. **Shared-versus-context decomposition:** stability of shared edges across bootstraps and reproducibility of context deviations.

### Required baselines

- cancer-stratified differential expression;
- covariate-adjusted linear model;
- PROGENy fixed weighted score;
- linear network diffusion;
- unconstrained MLP with matched parameter count;
- same architecture with degree-preserving randomised network;
- shared-only model without context gates;
- context-only model without shared propagation.

The neural model advances only if it improves leave-one-cancer-out effect recovery or independent phosphosite alignment, not merely within-cohort reconstruction error.

### Metrics

- correlation and signed concordance of predicted versus observed KRAS gene effects;
- precision/recall for held-out significant target direction;
- calibration of effect magnitude;
- CPTAC ERK-module association and site-level direction;
- edge-selection stability and pathway-distance enrichment;
- performance delta against linear diffusion and matched MLP.

## Immediate implementation gates

1. Reproduce TCGA gene-level effects with propensity weighting/matching, not regression alone.
2. Add major co-alterations and tumour purity to sensitivity analyses.
3. Test whether context-specific effects replicate between TCGA and CPTAC RNA within each cancer.
4. Construct the signed RAS/MAPK graph from frozen resources and quantify target coverage before model coding.
5. Implement linear diffusion and shared-only sparse neural propagation before context gating.

## Stop rules

- Stop if TCGA gene-level context effects do not replicate in CPTAC RNA beyond the shared PROGENy module.
- Stop if a linear covariate model or linear diffusion matches the context-gated network on transfer and phosphosite alignment.
- Stop if performance depends on including evaluated output genes in the context encoder.
- Stop any causal wording if external perturbation validation is absent or contradictory.

## Current status

The empirical feasibility gate is passed for a cross-cancer KRAS context project: RNA support exists in three cancers and independent phosphosite support in two. Model development is authorised only after the TCGA-to-CPTAC RNA replication and graph-coverage gates above.
