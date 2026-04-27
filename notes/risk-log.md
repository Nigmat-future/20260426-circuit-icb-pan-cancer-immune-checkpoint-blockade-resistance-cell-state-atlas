# Risk Log

| Date | Risk | Severity | Mitigation | Owner |
| --- | --- | --- | --- | --- |
| 2026-04-26 | Patient response labels differ across cohorts | High | Harmonize to controlled vocabulary and keep original labels | Analysis |
| 2026-04-26 | Cell-level tests can inflate significance | High | All response tests aggregate to patient/sample level | Analysis |
| 2026-04-26 | Bulk state score may reflect purity or immune infiltration | High | Adjust for purity proxies and benchmark against immune signatures | Analysis |
| 2026-04-26 | Large RDS files exceed local storage or time budget | Medium | Stage download; start with metadata, pseudobulk, and selected objects | Engineering |
| 2026-04-26 | Drug reversal evidence can be overinterpreted | Medium | Frame as prioritization, require multi-evidence support | Manuscript |
| 2026-04-26 | TCGA is not ICB-treated | Medium | Use TCGA only for biological characterization and prognosis | Manuscript |
| 2026-04-26 | Cross-cancer model may be dominated by melanoma | Medium | Leave-one-cancer-type-out validation and per-cancer reporting | Analysis |
| 2026-04-26 | Licensing or access varies by source | Medium | Keep access log and exclude controlled data from MVP unless approved | Governance |
