# ARISE execution correction: return to alteration-resolved attribution

## Status of the previous result

The previous global WP4155 expression-reconstruction prototype is an off-protocol auxiliary experiment. Its failure is valid for that architecture and task, but it is **not** a no-go decision for ARISE.

It deviated from the frozen proposal in three ways:

1. it evaluated pooled expression reconstruction rather than source-macro alteration-removal attribution;
2. it used pathway-module-to-PROGENy footprint edges as the primary decoder rather than signed TF-target output edges;
3. it did not evaluate each alteration through its own reachable nodes, carrier contrast, uncertainty and abstention rule.

The corresponding no-go language in `GO_NO_GO_2026-10-01.md` is superseded by this correction. Those results remain as a documented negative control and must not be deleted.

## Restored primary task

For each eligible alteration source `a` and patient `i`:

`phi(i,a) = f(A_i, z_i; G) - f(A_i without a, z_i; G)`

The model and evaluation must preserve:

- the alteration's functional direction at its source node;
- every signed directed propagation step;
- explicit intermediate-node states;
- signed TF-target decoding;
- source-reachable trans-expression only;
- source-specific carrier/non-carrier support;
- grouped attribution or abstention when the source is not identifiable.

## Execution order

### Gate A — source-specific graph manifest

For each source, freeze:

- functional alteration definition;
- source injection node and sign;
- reachable signalling nodes;
- reachable TFs;
- direct signed TF targets;
- paths with conflicting signs;
- evaluable protein/phosphosite intermediates.

Primary initial sources are KRAS canonical activating mutations and CTNNB1 activating mutations. PIK3R1 loss is withheld until its event-to-node sign is resolved; TP53 is retained as a secondary source because conservative CPTAC support is small.

### Gate B — observed direction audit

Within leakage-safe training data, estimate covariate-adjusted source-to-expression effects for every preregistered reachable direct target. Compare observed and graph-predicted signs separately for each source. PROGENy footprints are secondary summaries and cannot define direct target edges.

### Gate C — additive source-specific baselines

Fit covariates-only, direct alteration regression and linear signed propagation. Evaluation is source-macro and restricted to source-reachable trans targets. A pooled all-gene R2 is not a primary endpoint.

### Gate D — alteration-removal attribution

Only after Gates A–C pass, fit the smallest signed recurrent model. For each held-out carrier, toggle exactly one source and calculate `phi(i,a)`. Compare the carrier-average predicted contribution with held-out observed carrier contrasts, graph signs and locked intermediate evidence.

### Gate E — context admission

Context gates are admitted only if additive and linear signed models leave stable, preregistered non-additivity. No unconstrained patient embedding or RNA-derived gate is allowed.

## Stop rule

Stop the alteration-resolved project only if source-specific reachable targets lack incremental signal, real signed paths do not outperform matched path nulls, or alteration-removal attribution is unstable. Failure of pooled expression reconstruction alone is insufficient.
