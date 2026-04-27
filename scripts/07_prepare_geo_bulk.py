from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.geo_bulk import (
    prepare_gse115821,
    prepare_gse126044,
    prepare_gse135222,
    prepare_gse176307,
    prepare_gse78220,
    prepare_gse91061,
)
from circuit_icb.io import write_table


GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
HGNC_COMPLETE_SET_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"


def ensure_gene_info(path: Path) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(GENE_INFO_URL, path)


def ensure_hgnc_complete_set(path: Path) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(HGNC_COMPLETE_SET_URL, path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare locked GEO bulk ICB cohorts for CIRCUIT-ICB scoring.")
    parser.add_argument(
        "--cohort",
        required=True,
        choices=["GSE78220", "GSE91061", "GSE115821", "GSE176307", "GSE135222", "GSE126044"],
    )
    parser.add_argument("--out-dir", default="data/processed/bulk_icb")
    parser.add_argument("--gene-info", default="data/raw/ncbi/Homo_sapiens.gene_info.gz")
    parser.add_argument("--hgnc-complete-set", default="data/raw/hgnc/hgnc_complete_set.txt")
    args = parser.parse_args()

    out_dir = Path(args.out_dir) / args.cohort
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.cohort == "GSE78220":
        expression, clinical = prepare_gse78220(
            "data/raw/geo/GSE78220/GSE78220_PatientFPKM.xlsx",
            "data/raw/geo/GSE78220/GSE78220_series_matrix.txt.gz",
        )
    elif args.cohort == "GSE91061":
        ensure_gene_info(Path(args.gene_info))
        expression, clinical = prepare_gse91061(
            "data/raw/geo/GSE91061/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz",
            "data/raw/geo/GSE91061/GSE91061_series_matrix.txt.gz",
            args.gene_info,
            pretreatment_only=True,
        )
    elif args.cohort == "GSE115821":
        expression, clinical = prepare_gse115821(
            "data/raw/geo/GSE115821/GSE115821_MGH_counts.csv.gz",
            [
                "data/raw/geo/GSE115821/GSE115821-GPL11154_series_matrix.txt.gz",
                "data/raw/geo/GSE115821/GSE115821-GPL18573_series_matrix.txt.gz",
            ],
            pretreatment_only=True,
            collapse_to_patient=True,
        )
    elif args.cohort == "GSE176307":
        expression, clinical = prepare_gse176307(
            "data/raw/geo/GSE176307/GSE176307_salmon_tpm_gene.matrix.tsv.gz",
            "data/raw/geo/GSE176307/GSE176307_BACI_Omniseq_Sample_Name_Key_submitted_GEO_v2.csv.gz",
            "data/raw/geo/GSE176307/GSE176307_series_matrix.txt.gz",
        )
    elif args.cohort == "GSE135222":
        ensure_hgnc_complete_set(Path(args.hgnc_complete_set))
        expression, clinical = prepare_gse135222(
            "data/raw/geo/GSE135222/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",
            "data/raw/geo/GSE135222/GSE135222_series_matrix.txt.gz",
            args.hgnc_complete_set,
        )
    elif args.cohort == "GSE126044":
        expression, clinical = prepare_gse126044(
            "data/raw/geo/GSE126044/GSE126044_counts.txt.gz",
            "data/raw/geo/GSE126044/GSE126044_series_matrix.txt.gz",
        )
    else:
        raise ValueError(args.cohort)

    write_table(expression, out_dir / "expression.tsv")
    write_table(clinical, out_dir / "clinical.tsv")
    print(f"Wrote {args.cohort}: expression {expression.shape}, clinical {clinical.shape}")


if __name__ == "__main__":
    main()
