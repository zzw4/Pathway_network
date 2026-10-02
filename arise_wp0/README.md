# ARISE WP0

This directory contains the two-week feasibility audit for the frozen ARISE proposal.

The current stage is data and identifiability audit only. Neural-network development is explicitly out of scope until the WP0 gates pass.

## Reproducible first audit

The first script audits an existing local TCGA-BRCA cBioPortal-style study package without modifying it:

```powershell
python scripts/audit_tcga_legacy.py `
  --study-dir "C:\path\to\brca_tcga_pub2015" `
  --output-dir outputs/tcga_legacy
```

It reports:

- sample and patient counts by modality;
- strict patient intersections;
- clinical ER/PR/HER2 availability using IHC/FISH rather than PAM50;
- PIK3CA activating-hotspot and TP53 loss-of-function carrier counts;
- PIK3CA × TP53 2×2 support cells;
- ERBB2 and CCND1 high-level amplification counts;
- provenance and SHA-256 hashes for audited source files.

This legacy package is suitable for an early support audit, not yet the frozen training dataset. GDC-harmonized RNA/CNV/mutation data and CPTAC schema/overlap remain required.

