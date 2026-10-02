# ARISE multi-omics execution amendment

Date: 2 October 2026. Authorised by the user's instruction to continue with multi-omics.

## Scientific main line

Functional mutation/CNV changes a source state; typed signed relations propagate changes through signalling nodes and TFs; the model predicts downstream expression. Protein and phosphosite measurements constrain selected intermediate states during development. Alteration-removal contributions remain model-based explanations, not biological causal effects.

The cancer-specific network must cover the relevant WikiPathways collection, not just two successful mechanisms. WP4155 is an initial graph component; completeness of the cancer-specific collection is still unaudited. UCEC remains the existing feasibility workspace, not a claim that all cancers or datasets are settled.

## Changes to the historical plan

1. Development protein and phosphosite measurements may supervise explicit intermediate states.
2. Already inspected CPTAC-UCEC 2020 data are development/exploration data, not a locked external test.
3. Stable non-additivity is optional evidence, not the sole justification for neural modelling. Shared propagation parameters and intermediate supervision must be tested against non-neural counterparts.
4. Direction audits are descriptive until patient-level resampling/permutation validates them. Previous independent-gene binomial P values cannot support source pass/fail decisions.
5. Molecular subtype adjustment requires definition/provenance review; alteration-defined or RNA-derived labels are excluded from the primary adjustment set. Histology, purity and technical covariates require separate causal/predictive justification. No PAM50.
6. Intermediate and pathway parameters are not assumed identifiable from good RNA prediction. Interpretations are limited to evidence-supported nodes and stable contributions.

## Data roles and verified inventory

- TCGA-UCEC: mutation/CNV/RNA development; fold-specific transformations.
- CPTAC-UCEC 2020: paired multi-omics development. API assay lists contain 95 patients for mutation, CNA, RNA, protein and phosphoprotein. This is availability metadata, not a feature-level completeness count. The inconsistent `all` sample list has only 81 patients and must not define inclusion.
- CPTAC-UCEC 2023: candidate independent evaluation cohort. The NCI publication page reports 138 tumors and 20 normal tissues across ten platforms. Patient non-overlap, modality intersections and processed-matrix access remain to be verified. Do not inspect alteration-outcome associations while choosing this cohort.
- `uec_cptac_gdc`: API contains 442 samples from 241 patients and genomic/RNA profiles only. It is not automatically independent of either publication cohort. Tumor/normal filtering, duplicate resolution and case crosswalk are required.

Sources: https://gdc.cancer.gov/about-data/publications/CPTAC-3_2023_2 ; https://gdc.cancer.gov/about-data/publications/CPTAC-3_2020_4

## Node and observation semantics

One biological entity retains its identity. Only entities whose mechanisms require it receive distinct abundance and activity states; no blanket multiplication into one copy per omics layer.

Each edge records source state, target state, sign, mechanism, provenance and confidence. Kinase activation, transcriptional regulation, inhibitory phosphorylation, complex membership and protein degradation are not interchangeable edge types. Unresolved complex/neutral relations are excluded from signed propagation until represented explicitly. A database with no recorded conflict is not proof of biological consistency.

Examples:

- MEK activity -> ERK activity: activation.
- ELK1 activity -> FOS abundance: transcriptional activation.
- beta-catenin/TCF complex -> target RNA: complex-mediated regulation; CTNNB1 abundance is not automatically TCF activity.
- FOXO phosphorylation -> FOXO activity: sign requires site-specific annotation.

Observation heads map abundance states to total protein and activity states to curated activation/inhibition sites or substrate panels. A phosphosite without functional annotation is an observed feature, not an activity label. For sites affected by total abundance, a head models site intensity conditional on predicted protein abundance; protein adjustment is fitted inside training folds. No missing observation is encoded as zero.

## First model specification

Inputs at inference: eligible functional mutations, focal CNVs and admissible clinical/technical context. Protein/phosphoprotein are training supervision in the primary model, not clamped test inputs. Thus source removal can propagate through the complete learned chain; clamping a measured downstream mediator would obscure that contribution.

The shared small signed recurrent network produces target RNA and selected intermediate predictions. A nuisance branch receives covariates only. Each modality loss is averaged over observed entries per patient and then across patients, so large RNA panels do not dominate solely by feature count. All normalisation, observation filtering and loss weighting use development training folds only.

Loss: RNA prediction + weighted protein observation loss + weighted phosphosite observation loss + edge regularisation. Architecture/sign choices are audited before fitting. Observation heads do not establish that unobserved nodes are correct.

TCGA contributes RNA loss; paired CPTAC development patients contribute RNA and available intermediate losses. Cohort-specific measurement offsets are allowed but cannot receive alteration labels. Cross-cohort performance must be compared with cohort-aware non-neural controls.

## Baselines and claims

Compare covariates-only; direct ridge/reduced-rank alteration regression; linear signed propagation; linear multi-task propagation with the same intermediate supervision; equal-capacity unstructured neural network; signed neural network with and without intermediate supervision.

Graph controls include typed degree/sign/distance-aware rewiring and separately sign/edge ablations. Node labels and observation maps must be treated explicitly when constructing graph nulls.

Primary questions: does alteration information improve held-out reachable trans-RNA prediction, does intermediate supervision improve hidden intermediate prediction and stability, and does the real graph improve performance beyond matched controls? RNA accuracy alone cannot validate a path.

Evaluation uses patient-grouped nested CV. Held-out patients contribute no RNA/protein/phosphosite training values. A separate held-out-site task tests observation interpolation and is never called independent patient validation. Patient bootstrap estimates uncertainty; patient-label or residual permutations retain target correlations and covariate structure. Confounded source contrasts require overlap checks and sensitivity analyses rather than an unconditional 50% null.

## Immediate execution and completion criteria

1. Schema inventory (completed): `outputs/multiomics_schema_audit.json`, produced by `scripts/audit_multiomics_schema.py`.
2. Acquire processed development matrices and generate actual patient-by-modality and feature-missingness manifests; mutations absent from records count as wild type only for sequenced samples.
3. Establish publication-to-GDC-to-PDC crosswalk for the 2023 cohort and verify no development-patient overlap. Until complete, external evaluation feasibility is unresolved.
4. Build a mechanism/state audit for all core edges and extend from WP4155 to the cancer-specific pathway collection. Start with annotated observation heads; do not assign arbitrary sites to every node.
5. Repair source direction statistics and covariate provenance; retain KRAS and CTNNB1 as candidate cases pending corrected tests.
6. Fit CPU baselines and test actual effective sample size before granting a neural training budget. With roughly 95 paired development patients, a large free-edge model is unjustified; favour strongly shared weights and few observation-head parameters.

Next decision requires actual paired matrices, usable intermediate observations, credible independent validation access, and incremental held-out signal. No source swapping based on significant targets or external effect inspection.
