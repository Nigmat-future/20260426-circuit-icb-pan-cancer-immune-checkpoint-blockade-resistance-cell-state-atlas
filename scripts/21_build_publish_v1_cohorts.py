from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.benchmark_signatures import benchmark_summary, write_benchmark_signatures
from circuit_icb.cohort_schema import build_cohort_curation_row, standardize_publish_clinical, validate_publish_clinical
from circuit_icb.io import read_table, write_table


def standardize_known_cohort(cohort_id: str, expression_path: str, clinical_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = read_table(expression_path)
    clinical = read_table(clinical_path)
    if cohort_id == "GSE78220":
        clinical = standardize_publish_clinical(
            clinical,
            cohort_id="GSE78220",
            cancer_type="melanoma",
            therapy_class="anti-PD-1",
            sample_timing_col="treatment_timing",
            source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220",
            processing_note="Processed FPKM from GEO supplementary workbook; pre-treatment baseline samples.",
        )
    elif cohort_id == "GSE91061":
        clinical = standardize_publish_clinical(
            clinical,
            cohort_id="GSE91061",
            cancer_type="melanoma",
            therapy_class="anti-PD-1",
            sample_timing_col="treatment_timing",
            source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061",
            processing_note="Processed FPKM from GEO; pretreatment-only samples after Entrez-to-symbol mapping.",
        )
    elif cohort_id == "GSE115821":
        clinical = standardize_publish_clinical(
            clinical,
            cohort_id="GSE115821",
            cancer_type="melanoma",
            therapy_class="mixed anti-PD-1/anti-CTLA-4",
            sample_timing_col="treatment_timing",
            source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115821",
            processing_note="featureCounts matrix converted to log2(CPM+1); pretreatment samples collapsed to patient level.",
        )
    elif cohort_id == "GSE176307":
        clinical = standardize_publish_clinical(
            clinical,
            cohort_id="GSE176307",
            cancer_type="urothelial_carcinoma",
            therapy_class="anti-PD-1/anti-PD-L1",
            sample_timing_col="treatment_timing",
            source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176307",
            processing_note="Transparent GEO extension: Salmon gene TPM converted to log2(TPM+1); samples without matched GEO clinical metadata removed.",
        )
    elif cohort_id == "GSE135222":
        clinical = standardize_publish_clinical(
            clinical,
            cohort_id="GSE135222",
            cancer_type="non_small_cell_lung_cancer",
            therapy_class="anti-PD-1/anti-PD-L1",
            sample_timing_col="treatment_timing",
            source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222",
            processing_note="Transparent GEO extension: Ensembl IDs mapped through HGNC; expression converted to log2(value+1); DCB/NDB labels derived from PFS >= 180 days when RECIST response was unavailable.",
        )
    elif cohort_id == "GSE126044":
        clinical = standardize_publish_clinical(
            clinical,
            cohort_id="GSE126044",
            cancer_type="non_small_cell_lung_cancer",
            therapy_class="anti-PD-1",
            sample_timing_col="treatment_timing",
            source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044",
            processing_note="Transparent GEO extension: count matrix converted to log2(CPM+1); GEO pre-treatment patient_response labels used.",
        )
    else:
        raise ValueError(f"Unsupported cohort for publish_v1 standardization: {cohort_id}")
    return expression, clinical


def main() -> None:
    parser = argparse.ArgumentParser(description="Build publish_v1 cohort registry and standardized clinical tables.")
    parser.add_argument("--out-dir", default="data/processed/publish_v1")
    parser.add_argument("--benchmark-signatures", default="config/benchmark_signatures.yml")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_benchmark_signatures(args.benchmark_signatures)
    write_table(benchmark_summary(), "results/tables/publish_v1/benchmark_signature_summary.tsv")

    candidates = [
        ("GSE78220", "data/processed/bulk_icb/GSE78220/expression.tsv", "data/processed/bulk_icb/GSE78220/clinical.tsv"),
        ("GSE91061", "data/processed/bulk_icb/GSE91061/expression.tsv", "data/processed/bulk_icb/GSE91061/clinical.tsv"),
        ("GSE115821", "data/processed/bulk_icb/GSE115821/expression.tsv", "data/processed/bulk_icb/GSE115821/clinical.tsv"),
        ("GSE176307", "data/processed/bulk_icb/GSE176307/expression.tsv", "data/processed/bulk_icb/GSE176307/clinical.tsv"),
        ("GSE135222", "data/processed/bulk_icb/GSE135222/expression.tsv", "data/processed/bulk_icb/GSE135222/clinical.tsv"),
        ("GSE126044", "data/processed/bulk_icb/GSE126044/expression.tsv", "data/processed/bulk_icb/GSE126044/clinical.tsv"),
    ]
    registry_rows = []
    curation_rows = []
    summaries = []
    for cohort_id, expression_path, clinical_path in candidates:
        if not Path(expression_path).exists() or not Path(clinical_path).exists():
            curation_rows.append(
                build_cohort_curation_row(
                    cohort_id=cohort_id,
                    accession=cohort_id,
                    status="excluded",
                    reason="prepared expression or clinical table missing",
                    expression_file=expression_path,
                    clinical_file=clinical_path,
                    source_url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={cohort_id}",
                )
            )
            continue
        expression, clinical = standardize_known_cohort(cohort_id, expression_path, clinical_path)
        cohort_dir = out_dir / cohort_id
        write_table(expression, cohort_dir / "expression.tsv")
        write_table(clinical, cohort_dir / "clinical.tsv")
        summary = validate_publish_clinical(clinical)
        summaries.append(summary)
        registry_rows.append(
            {
                "cohort_id": cohort_id,
                "accession": cohort_id,
                "status": "included",
                "expression_file": str(cohort_dir / "expression.tsv"),
                "clinical_file": str(cohort_dir / "clinical.tsv"),
                "source_url": clinical["source_url"].iloc[0],
            }
        )
        curation_rows.append(
            build_cohort_curation_row(
                cohort_id=cohort_id,
                accession=cohort_id,
                status="included",
                reason="processed expression and response labels available",
                expression_file=str(cohort_dir / "expression.tsv"),
                clinical_file=str(cohort_dir / "clinical.tsv"),
                source_url=clinical["source_url"].iloc[0],
                notes=clinical["processing_note"].iloc[0],
            )
        )

    curation_rows.extend(
        [
            build_cohort_curation_row(
                cohort_id="IMvigor210",
                accession="IMvigor210",
                status="excluded_runtime_blocked",
                reason="official R package downloaded and inspected, but cds.RData is a Bioconductor S4 CountDataSet unsupported by pyreadr/rdata here and no Rscript runtime is available",
                source_url="http://research-pub.gene.com/IMvigor210CoreBiologies/",
                notes="Official tarball retained under data/raw/imvigor210; do not use non-transparent mirrors. Add when R/Bioconductor export is available.",
            ),
            build_cohort_curation_row(
                cohort_id="Gide2019_PRJEB23709",
                accession="PRJEB23709",
                status="excluded",
                reason="processed expression table not identified; raw FASTQ/SRA processing is outside the locked assumptions",
                source_url="https://www.ebi.ac.uk/ena/browser/view/PRJEB23709",
                notes="Can be added only if a transparent processed expression matrix and response metadata are supplied.",
            ),
        ]
    )

    registry = pd.DataFrame(registry_rows)
    curation = pd.DataFrame(curation_rows)
    summary = pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame()
    write_table(registry, "data/manifests/publish_v1_bulk_icb_cohorts.tsv")
    write_table(curation, "results/tables/publish_v1/cohort_curation_log.tsv")
    write_table(summary, "results/tables/publish_v1/cohort_summary.tsv")
    print(f"Wrote publish_v1 registry with {registry.shape[0]} included cohorts")


if __name__ == "__main__":
    main()
