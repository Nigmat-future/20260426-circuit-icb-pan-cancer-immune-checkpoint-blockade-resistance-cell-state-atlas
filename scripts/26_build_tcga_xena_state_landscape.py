from __future__ import annotations

import argparse
import gzip
import sys
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.geo_bulk import read_hgnc_ensembl_map
from circuit_icb.io import read_table, write_table
from circuit_icb.signature import load_signatures
from circuit_icb.tcga_xena import (
    cox_by_cancer_state,
    immune_proxy_correlations,
    summarize_state_landscape,
    zscore_signature_scores,
)


XENA_EXPRESSION_URL = "https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz"
XENA_HOST = "https://toil.xenahubs.net"
XENA_DATASET = "TcgaTargetGtex_rsem_gene_tpm"
XENA_PHENOTYPE_URL = "https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz"
XENA_SURVIVAL_URL = "https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp"
HGNC_COMPLETE_SET_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"


def ensure_download(url: str, path: str | Path) -> Path:
    out = Path(path)
    if out.exists() and out.stat().st_size > 0:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, out)
    return out


def load_tcga_metadata(phenotype_path: str | Path, survival_path: str | Path) -> pd.DataFrame:
    phenotype = pd.read_csv(phenotype_path, sep="\t", compression="infer", dtype=str, encoding="latin1")
    survival = pd.read_csv(survival_path, sep="\t", dtype=str, encoding="latin1")
    phenotype = phenotype.loc[
        phenotype["_study"].eq("TCGA") & phenotype["_sample_type"].eq("Primary Tumor")
    ].copy()
    merged = phenotype.merge(survival, on="sample", how="inner", suffixes=("", "_survival"))
    merged = merged.rename(columns={"cancer type abbreviation": "cancer_type"})
    return merged


def stream_signature_gene_expression(
    *,
    expression_url: str,
    hgnc_path: str | Path,
    samples: list[str],
    needed_symbols: set[str],
    out_path: str | Path,
) -> pd.DataFrame:
    out = Path(out_path)
    if out.exists() and out.stat().st_size > 0:
        return pd.read_csv(out, sep="\t", index_col=0, compression="infer")

    mapping = read_hgnc_ensembl_map(hgnc_path)
    ens_to_symbol = dict(zip(mapping["ensembl_gene_id"], mapping["symbol"].str.upper()))
    sample_set = set(samples)
    rows: list[list[object]] = []
    keep_samples: list[str] = []
    keep_positions: list[int] = []

    with urllib.request.urlopen(expression_url, timeout=600) as response:
        with gzip.GzipFile(fileobj=response) as gz:
            header = gz.readline().decode("utf-8").rstrip("\n").split("\t")
            for idx, sample in enumerate(header[1:], start=1):
                if sample in sample_set:
                    keep_samples.append(sample)
                    keep_positions.append(idx)
            try:
                for raw in gz:
                    parts = raw.decode("utf-8", errors="replace").rstrip("\n").split("\t")
                    ens = parts[0].split(".", 1)[0]
                    symbol = ens_to_symbol.get(ens)
                    if not symbol or symbol not in needed_symbols:
                        continue
                    values = [parts[pos] if pos < len(parts) else "" for pos in keep_positions]
                    rows.append([symbol, *values])
            except (EOFError, TimeoutError, OSError) as exc:
                print(f"Warning: Xena stream interrupted ({type(exc).__name__}); continuing with parsed signature genes if sufficient.", file=sys.stderr)

    if not rows:
        raise RuntimeError("No requested signature genes were found in the Xena expression matrix")
    expr = pd.DataFrame(rows, columns=["gene", *keep_samples])
    found_symbols = set(expr["gene"].astype(str).str.upper())
    if len(found_symbols) < max(10, int(len(needed_symbols) * 0.5)):
        raise RuntimeError(f"Only {len(found_symbols)}/{len(needed_symbols)} requested genes were parsed from Xena expression")
    expr[keep_samples] = expr[keep_samples].apply(pd.to_numeric, errors="coerce")
    expr = expr.groupby("gene", as_index=True)[keep_samples].mean()
    out.parent.mkdir(parents=True, exist_ok=True)
    expr.to_csv(out, sep="\t", compression="infer")
    return expr


def fetch_signature_gene_expression_api(
    *,
    samples: list[str],
    needed_symbols: set[str],
    out_path: str | Path,
) -> pd.DataFrame:
    out = Path(out_path)
    if out.exists() and out.stat().st_size > 0:
        return pd.read_csv(out, sep="\t", index_col=0, compression="infer")
    try:
        import xenaPython as xena
    except ImportError as exc:
        raise RuntimeError("xenaPython is required for API-based TCGA expression retrieval") from exc

    genes = sorted(needed_symbols)
    records = xena.dataset_gene_probe_avg(XENA_HOST, XENA_DATASET, samples, genes)
    rows = []
    for record in records:
        gene = str(record.get("gene", "")).upper()
        scores = record.get("scores", [])
        if not gene or not scores:
            continue
        values = scores[0] if isinstance(scores[0], list) else scores
        rows.append([gene, *values])
    if not rows:
        raise RuntimeError("No requested genes were returned by the Xena expression API")
    expr = pd.DataFrame(rows, columns=["gene", *samples])
    expr[samples] = expr[samples].apply(pd.to_numeric, errors="coerce")
    expr = expr.groupby("gene", as_index=True)[samples].mean()
    out.parent.mkdir(parents=True, exist_ok=True)
    expr.to_csv(out, sep="\t", compression="infer")
    return expr


def main() -> None:
    parser = argparse.ArgumentParser(description="Build TCGA Xena z-score state landscape for publish_v1.")
    parser.add_argument("--state-signatures", default="results/tables/state_annotation/final_state_signatures.yml")
    parser.add_argument("--benchmark-signatures", default="config/benchmark_signatures.yml")
    parser.add_argument("--raw-dir", default="data/raw/tcga_xena")
    parser.add_argument("--cache-dir", default="data/processed/tcga_xena")
    parser.add_argument("--out-dir", default="results/tables/publish_v1/tcga")
    args = parser.parse_args()

    raw_dir = Path(args.raw_dir)
    cache_dir = Path(args.cache_dir)
    out_dir = Path(args.out_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    phenotype = ensure_download(XENA_PHENOTYPE_URL, raw_dir / "TcgaTargetGTEX_phenotype.txt.gz")
    survival = ensure_download(XENA_SURVIVAL_URL, raw_dir / "Survival_SupplementalTable_S1_20171025_xena_sp.tsv")
    metadata = load_tcga_metadata(phenotype, survival)

    signatures = load_signatures(args.state_signatures)
    signatures.update(load_signatures(args.benchmark_signatures))
    needed_symbols = {gene for sig in signatures.values() for gene in (*sig.up, *sig.down)}
    samples = metadata["sample"].astype(str).tolist()
    try:
        expression = fetch_signature_gene_expression_api(
            samples=samples,
            needed_symbols=needed_symbols,
            out_path=cache_dir / "tcga_signature_gene_expression.tsv.gz",
        )
        expression_method = "xenaPython dataset_gene_probe_avg API"
    except Exception as exc:
        print(f"Warning: Xena API retrieval failed ({type(exc).__name__}: {exc}); falling back to streaming matrix.", file=sys.stderr)
        hgnc = ensure_download(HGNC_COMPLETE_SET_URL, Path("data/raw/hgnc/hgnc_complete_set.txt"))
        expression = stream_signature_gene_expression(
            expression_url=XENA_EXPRESSION_URL,
            hgnc_path=hgnc,
            samples=samples,
            needed_symbols=needed_symbols,
            out_path=cache_dir / "tcga_signature_gene_expression.tsv.gz",
        )
        expression_method = "streamed Xena gzip matrix with HGNC mapping"
    scores = zscore_signature_scores(expression, signatures, min_genes_present=3)
    scores = metadata.merge(scores, on="sample", how="inner")
    write_table(scores, out_dir / "tcga_state_scores.tsv")
    write_table(summarize_state_landscape(scores), out_dir / "tcga_state_landscape_by_cancer.tsv")
    write_table(cox_by_cancer_state(scores, time_col="OS.time", event_col="OS"), out_dir / "tcga_state_survival_os.tsv")
    write_table(cox_by_cancer_state(scores, time_col="PFI.time", event_col="PFI"), out_dir / "tcga_state_survival_pfi.tsv")
    write_table(immune_proxy_correlations(scores), out_dir / "tcga_state_immune_proxy_correlations.tsv")
    provenance = pd.DataFrame(
        [
            {"resource": "xena_expression_url", "value": XENA_EXPRESSION_URL},
            {"resource": "xena_host", "value": XENA_HOST},
            {"resource": "xena_dataset", "value": XENA_DATASET},
            {"resource": "expression_retrieval_method", "value": expression_method},
            {"resource": "xena_phenotype_url", "value": XENA_PHENOTYPE_URL},
            {"resource": "xena_survival_url", "value": XENA_SURVIVAL_URL},
            {"resource": "hgnc_complete_set_url", "value": HGNC_COMPLETE_SET_URL},
            {"resource": "n_tcga_primary_tumors_scored", "value": scores.shape[0]},
            {"resource": "n_signature_genes_extracted", "value": expression.shape[0]},
        ]
    )
    write_table(provenance, out_dir / "tcga_state_landscape_provenance.tsv")
    print(f"Wrote TCGA state landscape for {scores.shape[0]} primary tumors and {expression.shape[0]} signature genes")


if __name__ == "__main__":
    main()
