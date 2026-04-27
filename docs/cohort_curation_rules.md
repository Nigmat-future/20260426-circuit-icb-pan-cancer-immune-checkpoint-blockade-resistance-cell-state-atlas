# Cohort Curation Rules

## Required Cohort Fields

- cohort_id
- source accession or DOI
- cancer type
- treatment regimen
- treatment timing
- expression platform
- sample id
- patient id if available
- original response label
- harmonized response label
- survival endpoint and censoring if available
- covariates available

## Harmonized Response Rule

Primary binary endpoint:

- positive: CR, PR, responder, durable clinical benefit
- negative: SD, PD, non-responder, primary resistance
- excluded from primary binary analysis: mixed response, unknown

Sensitivity analyses:

- exclude on-treatment samples
- pretreatment-only models
- cancer-specific models

For cohorts without RECIST response, DCB/NDB may be derived from PFS only when the rule is explicit in `processing_note`; these cohorts must remain identifiable in the registry and should not be mixed silently with objective-response labels.

## Audit Trail

Keep the original label, harmonized label, and transformation rule in the same table. Never overwrite the source label.
