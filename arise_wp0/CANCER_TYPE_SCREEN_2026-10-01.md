# ARISE cancer-type screen — 1 October 2026

## Question

The ARISE method does not require breast cancer. We therefore screened the ten cancer types represented in the CPTAC pan-cancer resource rather than fixing BRCA in advance.

The first pass used TCGA PanCancer Atlas studies through the cBioPortal public REST API. It assessed conservative functional mutation classes, pair support and corresponding CPTAC patient counts. The second pass queried CPTAC molecular profiles and phosphoproteome metadata for candidate signalling intermediates.

## Ten-cancer mutation-support screen

Using broad nonsilent candidate-driver classes, the number of source pairs with all four cells at least 20 in TCGA was:

| Cancer | TCGA primary patients | Pairs passing 20-per-cell |
|---|---:|---:|
| UCEC | 529 | 13 |
| COADREAD | 594 | 12 |
| LUAD | 566 | 7 |
| HNSC | 523 | 7 |
| BRCA | 1,084 | 4 |
| GBM | 585 | 3 |
| CCRCC | 512 | 2 |
| LSCC | 487 | 0 |
| OV | 583 | 0 |
| PDAC | 184 | 0 |

Broad nonsilent counts are only a triage statistic. They were subsequently replaced by conservative functional classes for the four leading candidates.

## Functional-source refinement

### UCEC

Preferred pair: **PTEN loss-of-function + KRAS canonical hotspot**.

- TCGA cells (neither / KRAS only / PTEN only / both): 251 / 25 / 190 / 63.
- CPTAC functional carriers: PTEN 54, KRAS 29.
- CPTAC cells: 26 / 15 / 40 / 14.
- CPTAC linked mutation+CNA+RNA+protein+phosphoprotein samples: 95.
- PTEN→AKT/mTOR audit: all 11 candidate intermediate genes have phosphosite features; 192 candidate phosphosite features in total; all 11 proteins measured in at least 80% of samples.
- KRAS→RAF/MEK/ERK audit: all 8 candidate intermediate genes have phosphosite features; 83 candidate phosphosite features; all 8 proteins measured in at least 80% of samples.

Strength: two mechanistically distinct but potentially interacting signalling axes, adequate TCGA four-cell support, adequate CPTAC source carriers, and exceptionally complete intermediate coverage.

Critical risk: PTEN/KRAS status may be entangled with POLE, MSI, copy-number subtype, histology and grade. Subtype-adjusted local support is mandatory.

### LUAD

Preferred source candidates: **KRAS canonical hotspot + EGFR activating mutation**, with STK11 as a possible context variable rather than a primary source.

- TCGA KRAS/EGFR are essentially mutually exclusive, so their joint four-cell interaction is not identifiable.
- CPTAC carriers: KRAS 33, EGFR 32.
- CPTAC linked multi-omics samples: 110.
- KRAS route: 49 candidate phosphosite features across all 8 audited intermediates.
- EGFR route: 86 candidate phosphosite features across all 9 audited intermediates.

Strength: strongest and clearest receptor/RAS signalling biology with excellent external intermediate coverage.

Critical risk: mutual exclusivity makes this a source-comparison problem rather than an interaction/source-competition problem; smoking, ancestry and histology may dominate transportability.

### COADREAD

Preferred source candidates: **APC loss-of-function + KRAS canonical hotspot**.

- TCGA functional pair has adequate support.
- CPTAC carriers: APC 77, KRAS 27; CPTAC pair cells 28 / 5 / 55 / 22.
- CPTAC linked multi-omics samples: 95.
- APC/WNT route: 27 candidate phosphosite features, but only 2 of 6 audited proteins measured in at least 80% of samples.
- KRAS route: 43 candidate phosphosite features; 5 of 8 audited proteins measured in at least 80% of samples.

Strength: clear WNT and MAPK biology and strong TCGA support.

Critical risk: KRAS is strongly nested within APC-altered tumours, creating poor source-specific external support; APC-route total-protein coverage is weaker.

### BRCA

Preferred pair under conservative definitions: **PIK3CA canonical hotspot + TP53 damaging**.

- TCGA pair cells: 677 / 157 / 227 / 23.
- CPTAC carriers: PIK3CA 22, TP53 28; pair cells 79 / 21 / 15 / 7.
- CPTAC linked multi-omics samples: 122.
- PIK3CA route has excellent intermediate coverage.
- The audited TP53 route has only 7 phosphosite features across 4 measured genes and is not a comparably clear phosphorylation-propagation chain.

Strength: largest development and external cohorts among the four.

Critical risk: only one conservative TCGA pair passes the four-cell rule, CPTAC joint support is weak, and the two mechanisms have strongly unequal intermediate evidence.

## Provisional ranking

1. **UCEC — PIK3R1-LOF + KRAS-hotspot** (updated after covariate/RNA-footprint audit)
2. **LUAD — KRAS-hotspot + EGFR-activating**
3. **COADREAD — APC-LOF + KRAS-hotspot**
4. **BRCA — PIK3CA-hotspot + TP53-damaging**

The rankings answer different variants of the task:

- UCEC is best for source competition plus two distinct signalling axes.
- LUAD is best for clean pathway biology and molecular validation, but not for interaction between its two principal sources.
- COADREAD is biologically attractive but externally collinear.
- BRCA has the largest cohorts but the weakest second phospho mechanism under conservative source definitions.

## Decision

Breast cancer is no longer the default. **UCEC becomes the lead feasibility candidate**, with LUAD retained as the main alternative.

This is not yet a final cancer selection. UCEC advances only if candidate sources retain adequate propensity overlap and source-specific RNA increment after adjustment for molecular subtype, MSI/POLE status, histology, grade, purity and major co-alterations. CPTAC coverage counts must also be reduced to preregistered, directionally interpretable sites/modules with total-protein adjustment.

## UCEC confounding and RNA-footprint audit

The first covariate audit used 529 TCGA UCEC primary patients. Subtype was available for 507, grade and tumour type for all 529, and MANTIS MSI score for 526.

- PTEN truncating/splice loss is strongly associated with molecular subtype (Cramér's V = 0.55) and tumour type (V = 0.47). It must not be treated as an unconfounded source label.
- KRAS canonical hotspot status is also associated with subtype, but less strongly (V = 0.26).
- Within subtype, both sources remain represented in CN-low and MSI tumours, but the four PTEN/KRAS joint cells are too small for a subtype-stratified interaction claim. Therefore interaction is not a primary endpoint.

As a deliberately simple feasibility check, top-100 PROGENy MAPK and PI3K RNA footprints were calculated from TCGA expression and regressed on source status with subtype, tumour type, grade, MANTIS MSI score, and the alternate source as covariates (505 complete cases; HC3 robust errors).

- KRAS canonical hotspot retained the expected positive MAPK signal: beta = 1.42, 95% CI 0.79 to 2.04, P = 8.8e-6.
- PIK3R1 loss-of-function retained a positive PI3K signal: beta = 0.69, 95% CI 0.22 to 1.16, P = 0.0039.
- PTEN truncating/splice loss alone did not retain a convincing adjusted PI3K signal: beta = 0.33, P = 0.17.
- Adding PTEN deep deletion improved the broad functional definition only to beta = 0.44, P = 0.063; deep-deletion and biallelic-proxy groups each had only about 20 complete cases.
- PIK3CA canonical hotspots did not show the expected adjusted PI3K RNA footprint in this audit (P = 0.45).

This changes the leading UCEC source pair from PTEN/KRAS to **PIK3R1 loss-of-function / KRAS canonical hotspot**. TCGA four-cell counts are 373 / 68 / 68 / 20, while CPTAC UCEC counts are 53 / 25 / 13 / 4. Thus TCGA supports main effects and an exploratory interaction; CPTAC supports only separate pathway validation, not a four-cell interaction replication.

These PROGENy results are sanity checks, not the proposed model and not independent validation. Their value is to show that two mechanistically distinct alteration classes retain the expected expression footprints after major UCEC covariates are controlled. The next gate is whether curated, directionally interpretable CPTAC phosphosites show the matching intermediate effects after total-protein adjustment.

## CPTAC UCEC directional phosphosite gate

A small site panel was frozen before association testing: 12 PI3K/AKT/mTOR sites and four ERK activation-loop sites. Each phosphosite was adjusted for the corresponding total protein, genomic subtype, tumour purity, and the alternate source. The module score aggregates standardized protein-adjusted site residuals.

The proposed PIK3R1/KRAS dual-source design **failed the prespecified two-arm gate**.

- KRAS canonical hotspots passed the MAPK module check: beta = 0.459, 95% CI 0.082 to 0.836, two-sided P = 0.017. MAPK1 T185 and Y187 were individually positive after total-protein adjustment (P = 0.036 and 0.022, respectively).
- PIK3R1 loss-of-function failed the PI3K module check: beta = 0.001, 95% CI -0.291 to 0.293, P = 0.994.
- A bounded exploratory replacement screen using the same frozen PI3K module also failed for PIK3CA canonical hotspots (beta = -0.034, P = 0.842) and PTEN truncating/splice loss (beta = 0.018, P = 0.902).

Therefore UCEC should **not** proceed as a two-pathway PIK3R1/KRAS neural-network project. The KRAS-to-MAPK chain is supported independently at RNA and phosphoprotein levels, but the proposed PI3K arm is not. A defensible next design must either (a) become a cross-cancer, context-conditioned KRAS propagation project, or (b) identify a different cancer/source chain that passes the same RNA-plus-phosphosite gates. It must not rescue UCEC by post-hoc site selection.

## Generated machine-readable outputs

- `outputs/cancer_screen_tcga.json`
- `outputs/functional_source_screen.json`
- `outputs/functional_source_screen_cptac.json`
- `outputs/cptac_intermediate_coverage.json`
- `outputs/ucec_confounding_audit.json`
- `outputs/ucec_expression_signal.json`
- `outputs/ucec_cptac_phosphosite_gate.json`
