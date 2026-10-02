# Paired development baseline and entity graph repair

Controlling plan reread: `MULTIOMICS_EXECUTION_PLAN_2026-10-02.md`.

## Entity representation

`outputs/wp4155_restored_entity_graph.json` restores the original WikiPathways entity graph: 80 entities, including 3 complexes and 11 families/groups, with 28 original retained relations. The old 109 expanded gene edges and these 28 entity relations are different representations, not different counts of confirmed biochemical mechanisms.

Membership is stored separately, unsigned, and explicitly excluded from propagation. Family/complex membership can be traversed only for an uncertainty-labelled structural audit. No complex member-to-member causal edges are invented. Remaining relations retain unresolved state/mechanism labels pending curation; this is a representation repair, not a training-ready biological graph. Original data are preserved.

## Functional source feature package

Remote `outputs/paired_development/functional_source_features.npz` contains sample-aligned features and purity. Carriers: KRAS canonical 29; CTNNB1 exon-3 single-amino-acid substitutions at D32/S33/G34/S37/T41/S45 20; PIK3CA canonical 17; PTEN truncating/splice class 54; PIK3R1 truncating/splice class 17; TP53 truncating/splice class 7; TP53 canonical hotspot missense 5. Mutation classes preserve the historical class rules, whose functional confidence still needs review (particularly generic splice and Nonstop calls).

CTNNB1 encoding now requires a complete single-amino-acid substitution string, rather than a prefix match. The earlier count of 23 used broader matching; these counts must not be interchanged. Broader events remain unclassified, not automatically activating. PIK3R1 feature eligibility does not assert the direction of its node injection. CNV values remain available, but focality/cis exclusion is unresolved and no CNV pathway attribution is fitted here.

## Nested patient CV pilot

95 development patients, fixed five outer folds and three inner folds. Ridge alpha selected from 1/10/100 inside each outer training set. Feature imputation and input/output scaling are training-fold-specific. Protein/phosphosite measurements are outcomes, not predictors. Clinical control is cancer purity only; it does not remove all confounding.

| Endpoint | Full alteration model R2 | Corresponding-source increment R2 | Conditional patient-bootstrap 95% interval |
|---|---:|---:|---|
| KRAS-reachable direct RNA targets | 0.0540 | 0.0008 | -0.0031 to 0.0048 |
| CTNNB1-reachable direct RNA targets | 0.0612 | 0.0287 | 0.0175 to 0.0414 |
| CTNNB1 total protein | 0.4858 | 0.2476 | 0.0989 to 0.4118 |

The corresponding-source comparison refits and tunes a model with that feature omitted; it tests predictive information, not mutation-removal causal effects or a pathway-mediated contribution. Co-alterations stay in both models.

ERK phosphosite outcomes: full-model R2 -0.0004, increment beyond purity -0.0065. FOXO3 phosphosite outcomes: full-model R2 -0.0100, increment -0.0091. Neither has a stable gain in this small pilot. These raw site outcomes have not yet been protein-adjusted; they do not settle activity prediction feasibility.

## Limits and interpretation

This CPU pilot supports incremental CTNNB1 information at both protein and RNA endpoints. It does not establish that one mediates the other. KRAS-specific incremental RNA information is not established here. The paired cohort has only 95 patients; absent neural gain is not yet tested.

Target panels use the earlier component graph, pending semantic repair. This pilot uses features complete across the development cohort: 123/124 KRAS targets, 400/403 CTNNB1 targets, only two of four ERK sites and one of two FOXO3 sites. This global coverage screen makes it exploratory and unsuitable as the final leakage-safe benchmark; the final evaluator must filter inside training folds and score observed held-out entries with masks. No globally complete-case screen will be used for final model selection.

Bootstrap intervals condition on fitted fold predictions and fixed denominators; they do not include split/selection uncertainty. They are preliminary uncertainty summaries, not formal hypothesis tests. Broader covariate, missingness-aware and repeated-split sensitivity checks remain necessary before model admission.

Next: mechanism curation, source/covariate provenance, missingness-aware evaluation, independent cohort crosswalk, and the complete cancer-specific pathway inventory. Baseline results do not justify immediately training a large neural network.
