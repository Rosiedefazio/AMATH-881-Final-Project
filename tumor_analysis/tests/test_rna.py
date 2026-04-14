from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np
from scipy import io as scipy_io
from scipy import sparse

from tumor_atlas_ode.cli import main


def _write_gzip_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def _write_10x_sample(sample_dir: Path) -> None:
    genes = [
        "PHOX2B",
        "TH",
        "CHGA",
        "MKI67",
        "CD3D",
        "CD3E",
        "TRAC",
        "CD8A",
        "NKG7",
        "PRF1",
        "GZMB",
        "IFNG",
        "CCL5",
        "PDCD1",
        "LYZ",
        "FCN1",
        "S100A8",
        "S100A9",
        "CCR2",
        "B2M",
        "HLA-A",
        "TAP1",
        "CD68",
        "APOE",
        "C1QA",
        "MRC1",
        "IL10",
        "TGFB1",
    ]
    barcodes = ["cell_cancer", "cell_t", "cell_monocyte", "cell_macrophage"]
    matrix = np.zeros((len(genes), len(barcodes)), dtype=np.int32)

    def set_expr(cell_idx: int, gene_names: list[str], value: int = 25) -> None:
        for gene_name in gene_names:
            matrix[genes.index(gene_name), cell_idx] = value

    set_expr(0, ["PHOX2B", "TH", "CHGA", "MKI67"], value=40)
    set_expr(1, ["CD3D", "CD3E", "TRAC", "CD8A", "NKG7", "PRF1", "GZMB", "IFNG", "CCL5"], value=35)
    set_expr(1, ["PDCD1"], value=10)
    set_expr(2, ["LYZ", "FCN1", "S100A8", "S100A9", "CCR2", "B2M", "HLA-A", "TAP1"], value=30)
    set_expr(3, ["CD68", "APOE", "C1QA", "MRC1", "IL10", "TGFB1"], value=32)

    sample_dir.mkdir(parents=True, exist_ok=True)
    mm_buffer_path = sample_dir / "matrix.mtx"
    scipy_io.mmwrite(mm_buffer_path, sparse.csr_matrix(matrix))
    with mm_buffer_path.open("rb") as source, gzip.open(sample_dir / "matrix.mtx.gz", "wb") as target:
        target.write(source.read())
    mm_buffer_path.unlink()

    _write_gzip_text(
        sample_dir / "features.tsv.gz",
        "".join(f"gene_{idx}\t{name}\tGene Expression\n" for idx, name in enumerate(genes, start=1)),
    )
    _write_gzip_text(sample_dir / "barcodes.tsv.gz", "".join(f"{barcode}\n" for barcode in barcodes))


def test_summarize_scrna_cli_end_to_end(tmp_path: Path) -> None:
    publication_id = "demo-publication"
    sample_dir = (
        tmp_path
        / "data"
        / "raw"
        / publication_id
        / "downloads"
        / "scrna-seq"
        / "level-3"
        / "initial-diagnosis"
        / "participant-1"
        / "sample-1"
    )
    _write_10x_sample(sample_dir)

    exit_code = main(
        [
            "summarize-scrna",
            "--publication-id",
            publication_id,
            "--root",
            str(tmp_path),
        ]
    )

    assert exit_code == 0

    output_csv = tmp_path / "data" / "processed" / publication_id / "rna" / "sample_summaries.csv"
    assert output_csv.exists()

    with output_csv.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 1
    row = rows[0]
    assert row["timepoint_label"] == "initial-diagnosis"
    assert int(float(row["total_cells"])) == 4
    assert float(row["frac_cancer"]) > 0.0
    assert float(row["frac_cytotoxic_t"]) > 0.0
    assert float(row["frac_monocyte"]) > 0.0
    assert float(row["frac_macrophage"]) > 0.0
    assert float(row["interaction_mhci_cross_dressing"]) > 0.0
