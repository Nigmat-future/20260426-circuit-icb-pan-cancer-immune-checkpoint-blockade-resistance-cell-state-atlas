# Governance And Ethics

## Scope

This project uses public, de-identified research datasets. It is a retrospective computational study for biomarker and hypothesis generation. It does not provide patient-specific treatment recommendations.

## Data Access

All data sources must be logged with:

- URL or accession
- license or data-use statement
- date accessed
- downloaded file names
- checksums when available
- transformation steps

Controlled-access resources are excluded from the minimum publishable unit unless access is independently approved.

## Reproducibility

The project should release:

- manifests
- source code
- environment files
- workflow definitions
- exact parameter configuration
- release tables
- figure-generating notebooks or scripts

Large raw datasets are not committed to version control. The workflow stores download instructions and checksums.

## Privacy

No attempt should be made to re-identify patients. Clinical metadata should be harmonized only at the level necessary for response, survival, cancer type, treatment, and covariate analysis.

## Reporting

The manuscript should follow:

- transparent cohort inclusion/exclusion reporting
- TRIPOD-AI style reporting for prediction components where relevant
- FAIR principles for reusable resources
- open-source release expectations for computational resource papers

## Claim Boundaries

Allowed claims:

- malignant-cell state association with ICB response
- retrospective prediction or stratification performance
- prioritized candidate vulnerabilities
- perturbational hypotheses for combination therapy

Avoided claims:

- proven causal resistance mechanism
- clinically actionable treatment assignment
- validated therapeutic efficacy
- pan-cancer universality without leave-cancer-out support
