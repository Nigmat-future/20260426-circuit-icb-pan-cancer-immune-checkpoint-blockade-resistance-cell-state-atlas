# Target Scoring

## Evidence Components

Each target receives independent evidence values from 0 to 1:

- state_specificity: malignant resistance-state specificity
- clinical_association: higher in ICB non-response or poor outcome
- protein_evidence: CPTAC protein or phosphoprotein support
- dependency: DepMap dependency in relevant lineage or high-state models
- drug_sensitivity: PRISM/GDSC sensitivity association
- perturbation_reversal: LINCS/CLUE reversal of state signature
- safety: lower normal-tissue risk from GTEx/HPA
- druggability: availability of drug, inhibitor class, or tractable target family

## Default Formula

`target_priority_score = weighted_sum(evidence) / sum(weights)`

`conservative_priority_score = target_priority_score * evidence_completeness`

The conservative score penalizes targets with many missing evidence layers.

## Interpretation

Top targets are nominated as hypotheses. They are not validated therapies and should not be used for patient treatment decisions.
