# Paired development data and graph state audit

Controlling plan: `MULTIOMICS_EXECUTION_PLAN_2026-10-02.md`, reread before this work.

## Completed acquisition

Remote package: `/mnt/primary/arise_wp0/outputs/paired_development/paired_development.npz`.

The package has one sample per patient for 95 CPTAC-UCEC 2020 development patients. All 95 are in the sequenced list and have at least one finite observation in CNA, RNA, protein and the selected phosphosite panel. Missing mutation records are interpretable only within this sequenced and gene-profile scope.

| Matrix | Dimensions | Observed entries |
|---|---|---|
| RNA: core plus direct TF targets | 95 × 797 | 99.50% |
| Total protein: WP4155 core | 95 × 63 | 81.70% |
| GISTIC CNA: WP4155 core | 95 × 63 | 100% |
| Selected phosphosite panel | 95 × 20 | heterogeneous; see feature report |

The 760 direct-target genes cover all eligible terminal TF decoder edges in the current component, rather than only KRAS/CTNNB1 targets. The package preserves raw API values and missingness, contains no fitted scaling or imputation, and is a development resource. GISTIC gene values alone do not establish focality: focal CNV attribution still requires segment/region annotation. Mutation source functional classes are not yet encoded in the NPZ; detailed records are cached separately. This is a scoped component dataset, not a completed full cancer network dataset.

Dataset SHA256: `28bd9896b9d5076eeec52fd49a7be9a66784baedcf7e2fa780b3605a24e0da9b`.

## Intermediate observation feasibility

ERK has useful activation-site coverage: MAPK1 T185/Y187 in 61/95 patients and MAPK3 T202/Y204 in 55/95. Single-site values are not proof of dual-site phosphorylation. FOXO3 inhibitory T32/S253 sites cover 95/89 patients, with total protein in all 95. These are candidate observation heads subject to localisation and measurement-quality review.

CTNNB1 total protein covers all 95, but S33/S37/T41/S45 cover only 15/0/7/15. It supports an abundance observation head, not a well-covered phosphosite activity head. Alterations at these residues can affect peptide detection and must be checked before interpreting missingness or intensity. MYC has no protein observation in this scoped profile, despite being a major RNA decoder TF.

## Graph construction audit

Among 109 current core edges, automated audit flags overlap:

- 91 lack sufficiently resolved activity/abundance/mechanism semantics.
- 30 arise from interactions touching WikiPathways complex objects expanded into member-to-member edges.
- 3 are neutral or undirected.
- 4 CTNNB1-to-TCF/LEF relations require explicit complex representation.
- 2 core pairs have conflicting signs somewhere in candidate integrated evidence: AKT1 -> FOXO3 and TP53 -> POLK.

These are audit flags, not evidence that every flagged edge is false. Conversely, the earlier zero-conflict path result examined only retained paths and cannot establish that integrated biological evidence is conflict-free. No edge is automatically promoted to an activity mechanism merely from an arrow symbol.

## Next implementation dependencies

1. Curate complex nodes and edge states before signed neural propagation.
2. Verify phosphosite localisation/aggregation and assay scales; API values are convenient but not a substitute for original processed-data documentation.
3. Establish 2023 independent-cohort case crosswalk and modality intersection.
4. Construct functional source features, admissible covariates and patient splits; fit CPU baselines with fold-specific transforms.

No predictive, causal or neural-admission conclusion follows from this coverage audit.
