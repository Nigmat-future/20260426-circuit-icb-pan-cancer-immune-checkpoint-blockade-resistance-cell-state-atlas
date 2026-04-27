from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd


def _split_geo_line(line: str) -> list[str]:
    parts = line.rstrip("\n").split("\t")
    return [part.strip().strip('"') for part in parts[1:]]


def parse_series_matrix_metadata(path: str | Path) -> pd.DataFrame:
    """Parse sample metadata lines from a GEO series matrix file."""

    opener = gzip.open if str(path).endswith(".gz") else open
    titles: list[str] = []
    accessions: list[str] = []
    source_names: list[str] = []
    characteristics: list[list[str]] = []

    with opener(path, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_title"):
                titles = _split_geo_line(line)
            elif line.startswith("!Sample_geo_accession"):
                accessions = _split_geo_line(line)
            elif line.startswith("!Sample_source_name_ch1"):
                source_names = _split_geo_line(line)
            elif line.startswith("!Sample_characteristics_ch1"):
                characteristics.append(_split_geo_line(line))

    if not titles:
        raise ValueError(f"No !Sample_title line found in {path}")

    rows = []
    for idx, title in enumerate(titles):
        row = {
            "sample_title": title,
            "geo_accession": accessions[idx] if idx < len(accessions) else "",
            "source_name": source_names[idx] if idx < len(source_names) else "",
        }
        for characteristic in characteristics:
            if idx >= len(characteristic):
                continue
            value = characteristic[idx]
            if ":" not in value:
                continue
            key, raw = value.split(":", 1)
            key = re.sub(r"[^a-z0-9]+", "_", key.strip().lower()).strip("_")
            if key:
                row[key] = raw.strip()
        rows.append(row)
    return pd.DataFrame(rows)


def read_ncbi_gene_info(path: str | Path) -> pd.DataFrame:
    """Read NCBI Homo_sapiens.gene_info.gz as Entrez to symbol map."""

    columns = [
        "tax_id",
        "GeneID",
        "Symbol",
        "LocusTag",
        "Synonyms",
        "dbXrefs",
        "chromosome",
        "map_location",
        "description",
        "type_of_gene",
        "Symbol_from_nomenclature_authority",
        "Full_name_from_nomenclature_authority",
        "Nomenclature_status",
        "Other_designations",
        "Modification_date",
        "Feature_type",
    ]
    df = pd.read_csv(path, sep="\t", comment="#", names=columns, compression="infer", dtype={"GeneID": str})
    return df[["GeneID", "Symbol"]].drop_duplicates()


def _collapse_duplicate_genes(df: pd.DataFrame, gene_col: str) -> pd.DataFrame:
    values = df.dropna(subset=[gene_col]).copy()
    numeric_cols = [col for col in values.columns if col != gene_col]
    values[numeric_cols] = values[numeric_cols].apply(pd.to_numeric, errors="coerce")
    return values.groupby(gene_col, as_index=False)[numeric_cols].mean()


def prepare_gse78220(expression_xlsx: str | Path, series_matrix: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_excel(expression_xlsx, sheet_name="FPKM")
    expr = expr.rename(columns={"Gene": "gene"})
    expr = _collapse_duplicate_genes(expr, "gene")
    expression = expr.set_index("gene").T.reset_index().rename(columns={"index": "sample"})
    expression["sample_title"] = expression["sample"].str.replace(r"\.baseline$", "", regex=True)

    meta = parse_series_matrix_metadata(series_matrix)
    clinical = expression[["sample", "sample_title"]].merge(meta, on="sample_title", how="left")
    clinical["cohort_id"] = "GSE78220"
    clinical["response_original"] = clinical.get("anti_pd_1_response", "")
    clinical["treatment_timing"] = "pre-treatment"
    clinical["cancer_type"] = "melanoma"
    return expression.drop(columns=["sample_title"]), clinical


def prepare_gse91061(
    expression_csv_gz: str | Path,
    series_matrix: str | Path,
    gene_info_gz: str | Path,
    *,
    pretreatment_only: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(expression_csv_gz, compression="infer")
    first = expr.columns[0]
    expr = expr.rename(columns={first: "GeneID"})
    expr["GeneID"] = expr["GeneID"].astype(str)
    gene_info = read_ncbi_gene_info(gene_info_gz)
    expr = expr.merge(gene_info, on="GeneID", how="left")
    expr = expr.drop(columns=["GeneID"])
    expr = expr.rename(columns={"Symbol": "gene"})
    expr = _collapse_duplicate_genes(expr, "gene")
    expression = expr.set_index("gene").T.reset_index().rename(columns={"index": "sample"})

    meta = parse_series_matrix_metadata(series_matrix)
    clinical = expression[["sample"]].merge(meta, left_on="sample", right_on="sample_title", how="left")
    clinical["cohort_id"] = "GSE91061"
    clinical["response_original"] = clinical.get("response", "")
    clinical["treatment_timing"] = clinical.get("visit_pre_or_on_treatment", "")
    clinical["cancer_type"] = "melanoma"

    if pretreatment_only:
        keep = clinical["treatment_timing"].astype(str).str.lower().eq("pre")
        samples = clinical.loc[keep, "sample"].tolist()
        expression = expression.loc[expression["sample"].isin(samples)].copy()
        clinical = clinical.loc[keep].copy()
    return expression, clinical


def _clean_gse115821_sample(sample: str) -> str:
    sample = str(sample).strip()
    sample = re.sub(r"\.bam$", "", sample)
    return sample


def _log2_cpm_from_featurecounts(counts: pd.DataFrame) -> pd.DataFrame:
    sample_cols = [col for col in counts.columns if col not in {"Geneid", "Chr", "Start", "End", "Strand", "Length"}]
    values = counts[["Geneid", *sample_cols]].rename(columns={"Geneid": "gene"}).copy()
    values["gene"] = values["gene"].astype(str)
    values[sample_cols] = values[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    collapsed = values.groupby("gene", as_index=False)[sample_cols].sum()
    matrix = collapsed.set_index("gene")
    library_sizes = matrix.sum(axis=0).replace(0, np.nan)
    cpm = matrix.divide(library_sizes, axis=1) * 1_000_000.0
    log_cpm = np.log2(cpm + 1.0)
    return log_cpm.T.reset_index().rename(columns={"index": "sample"})


def _log2_cpm_from_gene_counts(counts: pd.DataFrame, gene_col: str) -> pd.DataFrame:
    sample_cols = [col for col in counts.columns if col != gene_col]
    values = counts[[gene_col, *sample_cols]].rename(columns={gene_col: "gene"}).copy()
    values["gene"] = values["gene"].astype(str)
    values[sample_cols] = values[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    collapsed = values.groupby("gene", as_index=False)[sample_cols].sum()
    matrix = collapsed.set_index("gene")
    library_sizes = matrix.sum(axis=0).replace(0, np.nan)
    cpm = matrix.divide(library_sizes, axis=1) * 1_000_000.0
    log_cpm = np.log2(cpm + 1.0)
    return log_cpm.T.reset_index().rename(columns={"index": "sample"})


def read_hgnc_ensembl_map(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    cols = {"ensembl_gene_id", "symbol"}
    if not cols.issubset(df.columns):
        raise ValueError(f"HGNC table must contain {cols}; found {df.columns.tolist()[:10]}")
    out = df[["ensembl_gene_id", "symbol"]].dropna().copy()
    out["ensembl_gene_id"] = out["ensembl_gene_id"].astype(str).str.split(".", n=1).str[0]
    out["symbol"] = out["symbol"].astype(str)
    return out.drop_duplicates("ensembl_gene_id")


def _dcb_from_pfs(pfs_time: object, pfs_event: object, *, threshold_days: int = 180) -> str:
    time = pd.to_numeric(pd.Series([pfs_time]), errors="coerce").iloc[0]
    event = pd.to_numeric(pd.Series([pfs_event]), errors="coerce").iloc[0]
    if pd.isna(time):
        return ""
    if time >= threshold_days:
        return "DCB"
    if not pd.isna(event) and int(event) == 1:
        return "NDB"
    return ""


def prepare_gse115821(
    counts_csv_gz: str | Path,
    series_matrices: list[str | Path],
    *,
    pretreatment_only: bool = True,
    collapse_to_patient: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(counts_csv_gz, compression="infer", encoding="utf-8-sig")
    expression = _log2_cpm_from_featurecounts(counts)
    expression["sample_title"] = expression["sample"].map(_clean_gse115821_sample)

    metadata = []
    for matrix in series_matrices:
        metadata.append(parse_series_matrix_metadata(matrix))
    meta = pd.concat(metadata, ignore_index=True)
    clinical = expression[["sample", "sample_title"]].merge(meta, on="sample_title", how="left")
    clinical["cohort_id"] = "GSE115821"
    clinical["response_original"] = clinical.get("response", "")
    clinical["treatment_timing"] = clinical.get("treatment_state", "")
    clinical["cancer_type"] = "melanoma"
    clinical["therapy_class"] = clinical.get("antibody", "")
    clinical["patient_id"] = clinical.get("patient_id", clinical["sample_title"]).astype(str).str.replace("MGH Patient ", "MGH", regex=False)
    clinical["timepoint_numeric"] = pd.to_numeric(clinical.get("timepoint", pd.NA), errors="coerce")

    if pretreatment_only:
        keep = clinical["treatment_timing"].astype(str).str.contains("PRE", case=False, na=False)
        clinical = clinical.loc[keep].copy()
        expression = expression.loc[expression["sample"].isin(clinical["sample"])].copy()

    expression = expression.drop(columns=["sample_title"])
    if collapse_to_patient and not clinical.empty:
        gene_cols = [col for col in expression.columns if col != "sample"]
        merged_expr = expression.merge(clinical[["sample", "patient_id"]], on="sample", how="inner")
        expression = merged_expr.groupby("patient_id", as_index=False)[gene_cols].mean().rename(columns={"patient_id": "sample"})
        clinical = (
            clinical.sort_values(["patient_id", "timepoint_numeric"], na_position="last")
            .drop_duplicates("patient_id")
            .copy()
        )
        clinical["sample"] = clinical["patient_id"]

    return expression, clinical


def prepare_gse176307(
    tpm_matrix_tsv_gz: str | Path,
    sample_key_csv_gz: str | Path,
    series_matrix: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(tpm_matrix_tsv_gz, sep="\t", compression="infer")
    first = expr.columns[0]
    expr = expr.rename(columns={first: "gene"})
    key = pd.read_csv(sample_key_csv_gz, compression="infer")
    key = key.rename(columns={"Sample ID": "sample", "Omniseq_RS_ID (RNAseq)": "rs_id"})
    rs_to_sample = dict(zip(key["rs_id"].astype(str), key["sample"].astype(str)))
    sample_cols = [col for col in expr.columns if col != "gene"]
    expr = expr.rename(columns={col: rs_to_sample.get(str(col), str(col)) for col in sample_cols})
    sample_cols = [col for col in expr.columns if col != "gene"]
    expr[sample_cols] = expr[sample_cols].apply(pd.to_numeric, errors="coerce")
    expr = _collapse_duplicate_genes(expr, "gene")
    expr[sample_cols] = np.log2(expr[sample_cols].clip(lower=0.0) + 1.0)
    expression = expr.set_index("gene").T.reset_index().rename(columns={"index": "sample"})

    meta = parse_series_matrix_metadata(series_matrix)
    meta["sample"] = meta["sample_title"].astype(str).str.extract(r"(BACI\d+)", expand=False)
    clinical = expression[["sample"]].merge(meta, on="sample", how="left")
    clinical["cohort_id"] = "GSE176307"
    clinical["response_original"] = clinical.get("io_response", "")
    clinical["treatment_timing"] = "pre-treatment"
    clinical["cancer_type"] = "urothelial_carcinoma"
    clinical["therapy_class"] = clinical.get("io_therapy", "immune checkpoint blockade")
    clinical["os_time"] = pd.to_numeric(clinical.get("overall_survival", pd.NA), errors="coerce")
    clinical["os_event"] = clinical.get("alive", "").astype(str).str.lower().map({"yes": 0, "no": 1})
    clinical["pfs_time"] = pd.to_numeric(clinical.get("pfs", pd.NA), errors="coerce")
    clinical["pfs_event"] = clinical.get("progressed", "").astype(str).str.lower().map({"yes": 1, "no": 0})
    clinical = clinical.loc[clinical["geo_accession"].notna()].copy()
    expression = expression.loc[expression["sample"].isin(clinical["sample"])].copy()
    return expression, clinical


def prepare_gse135222(
    expression_tsv_gz: str | Path,
    series_matrix: str | Path,
    hgnc_complete_set_tsv: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(expression_tsv_gz, sep="\t", compression="infer")
    expr = expr.rename(columns={expr.columns[0]: "ensembl_gene_id"})
    expr["ensembl_gene_id"] = expr["ensembl_gene_id"].astype(str).str.split(".", n=1).str[0]
    mapping = read_hgnc_ensembl_map(hgnc_complete_set_tsv)
    expr = expr.merge(mapping, on="ensembl_gene_id", how="left")
    expr = expr.drop(columns=["ensembl_gene_id"]).rename(columns={"symbol": "gene"})
    sample_cols = [col for col in expr.columns if col != "gene"]
    expr[sample_cols] = expr[sample_cols].apply(pd.to_numeric, errors="coerce")
    expr = _collapse_duplicate_genes(expr, "gene")
    sample_cols = [col for col in expr.columns if col != "gene"]
    expr[sample_cols] = np.log2(expr[sample_cols].clip(lower=0.0) + 1.0)
    expression = expr.set_index("gene").T.reset_index().rename(columns={"index": "sample"})

    meta = parse_series_matrix_metadata(series_matrix)
    meta["sample"] = meta["sample_title"].astype(str).str.replace(" ", "", regex=False)
    clinical = expression[["sample"]].merge(meta, on="sample", how="left")
    clinical["cohort_id"] = "GSE135222"
    clinical["treatment_timing"] = "pre-treatment"
    clinical["cancer_type"] = "non_small_cell_lung_cancer"
    clinical["therapy_class"] = "anti-PD-1/anti-PD-L1"
    clinical["pfs_time"] = pd.to_numeric(clinical.get("pfs_time", pd.NA), errors="coerce")
    clinical["pfs_event"] = pd.to_numeric(clinical.get("progression_free_survival_pfs", pd.NA), errors="coerce")
    clinical["response_original"] = [
        _dcb_from_pfs(time, event) for time, event in zip(clinical["pfs_time"], clinical["pfs_event"])
    ]
    return expression, clinical


def prepare_gse126044(counts_txt_gz: str | Path, series_matrix: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(counts_txt_gz, sep="\t", compression="infer")
    counts = counts.rename(columns={counts.columns[0]: "gene"})
    expression = _log2_cpm_from_gene_counts(counts, "gene")

    meta = parse_series_matrix_metadata(series_matrix)
    meta["sample"] = meta["sample_title"].astype(str).str.replace("RNA-seq_", "", regex=False)
    clinical = expression[["sample"]].merge(meta, on="sample", how="left")
    clinical["cohort_id"] = "GSE126044"
    clinical["response_original"] = clinical.get("patient_response", "")
    clinical["treatment_timing"] = clinical.get("treatment_status", "pre-treatment")
    clinical["cancer_type"] = "non_small_cell_lung_cancer"
    clinical["therapy_class"] = "anti-PD-1"
    return expression, clinical
