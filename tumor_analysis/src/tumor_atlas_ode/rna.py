from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy import io as scipy_io
from scipy import sparse

from tumor_atlas_ode.markers import CELL_TYPE_MARKERS, PROGRAM_MARKERS
from tumor_atlas_ode.utils import ensure_dir, write_csv


def discover_10x_sample_dirs(root: Path) -> list[Path]:
    sample_dirs: list[Path] = []
    for matrix_path in root.rglob("matrix.mtx.gz"):
        sample_dir = matrix_path.parent
        if (sample_dir / "features.tsv.gz").exists() and (sample_dir / "barcodes.tsv.gz").exists():
            sample_dirs.append(sample_dir)
    return sorted(sample_dirs)


def _read_tsv_column(path: Path, preferred_index: int = 1) -> list[str]:
    values: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            if preferred_index < len(parts):
                values.append(parts[preferred_index])
            else:
                values.append(parts[0])
    return values


def load_10x_matrix(sample_dir: Path) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    with gzip.open(sample_dir / "matrix.mtx.gz", "rt", encoding="utf-8") as handle:
        matrix = scipy_io.mmread(handle).tocsr().astype(np.float32)
    features = _read_tsv_column(sample_dir / "features.tsv.gz", preferred_index=1)
    barcodes = _read_tsv_column(sample_dir / "barcodes.tsv.gz", preferred_index=0)
    if matrix.shape[0] == len(barcodes) and matrix.shape[1] == len(features):
        matrix = matrix.transpose().tocsr()
    if matrix.shape[0] != len(features):
        raise RuntimeError(
            f"Feature count mismatch in {sample_dir}: matrix has {matrix.shape[0]} rows, features file has {len(features)} rows."
        )
    if matrix.shape[1] != len(barcodes):
        raise RuntimeError(
            f"Barcode count mismatch in {sample_dir}: matrix has {matrix.shape[1]} columns, barcodes file has {len(barcodes)} rows."
        )
    return matrix, features, barcodes


def log_normalize(matrix: sparse.csr_matrix, scale_factor: float = 1e4) -> sparse.csr_matrix:
    library_sizes = np.asarray(matrix.sum(axis=0)).ravel()
    safe_sizes = np.where(library_sizes > 0, library_sizes, 1.0)
    scale = scale_factor / safe_sizes
    normalized = matrix.multiply(scale)
    normalized.data = np.log1p(normalized.data)
    return normalized.tocsr()


def build_gene_index(features: Iterable[str]) -> dict[str, list[int]]:
    gene_index: dict[str, list[int]] = {}
    for idx, gene in enumerate(features):
        gene_index.setdefault(gene.upper(), []).append(idx)
    return gene_index


def panel_score(matrix: sparse.csr_matrix, gene_index: dict[str, list[int]], markers: list[str]) -> tuple[np.ndarray, int]:
    indices: list[int] = []
    for marker in markers:
        indices.extend(gene_index.get(marker.upper(), []))
    if not indices:
        return np.zeros(matrix.shape[1], dtype=np.float32), 0
    submatrix = matrix[indices, :]
    return np.asarray(submatrix.mean(axis=0)).ravel().astype(np.float32), len(indices)


def zscore(values: np.ndarray) -> np.ndarray:
    mean = float(values.mean())
    std = float(values.std())
    if std == 0.0:
        return np.zeros_like(values)
    return (values - mean) / std


def classify_cells(cell_type_scores: dict[str, np.ndarray]) -> np.ndarray:
    labels = list(cell_type_scores)
    stacked = np.vstack([zscore(cell_type_scores[label]) for label in labels])
    best_index = stacked.argmax(axis=0)
    best_score = stacked.max(axis=0)
    assigned = np.array([labels[idx] for idx in best_index], dtype=object)
    assigned[best_score <= 0.0] = "unknown"
    return assigned


def _mean_over_mask(values: np.ndarray, mask: np.ndarray) -> float:
    if not np.any(mask):
        return 0.0
    return float(values[mask].mean())


def parse_metadata_from_sample_dir(sample_dir: Path, downloads_root: Path) -> dict[str, str]:
    relative = sample_dir.relative_to(downloads_root)
    parts = relative.parts
    metadata = {
        "assay": parts[0] if len(parts) > 0 else "",
        "level": parts[1] if len(parts) > 1 else "",
        "timepoint_label": parts[2] if len(parts) > 2 else "",
        "participant_id": parts[3] if len(parts) > 3 else "",
        "sample_key": parts[4] if len(parts) > 4 else sample_dir.name,
        "sample_dir": str(sample_dir),
    }
    return metadata


def summarize_10x_sample(sample_dir: Path, downloads_root: Path) -> dict[str, float | str]:
    raw_matrix, features, _barcodes = load_10x_matrix(sample_dir)
    matrix = log_normalize(raw_matrix)
    gene_index = build_gene_index(features)

    cell_type_scores: dict[str, np.ndarray] = {}
    matched_markers: dict[str, int] = {}
    for label, markers in CELL_TYPE_MARKERS.items():
        score, matched = panel_score(matrix, gene_index, markers)
        cell_type_scores[label] = score
        matched_markers[label] = matched

    assignments = classify_cells(cell_type_scores)
    metadata = parse_metadata_from_sample_dir(sample_dir, downloads_root)
    summary: dict[str, float | str] = dict(metadata)
    summary["total_cells"] = int(raw_matrix.shape[1])

    class_counts = {label: int(np.sum(assignments == label)) for label in list(CELL_TYPE_MARKERS) + ["unknown"]}
    classified_total = max(sum(class_counts[label] for label in CELL_TYPE_MARKERS), 1)
    summary["classified_cells"] = classified_total

    for label, count in class_counts.items():
        summary[f"cells_{label}"] = count
    for label in CELL_TYPE_MARKERS:
        summary[f"frac_{label}"] = class_counts[label] / classified_total
        summary[f"matched_markers_{label}"] = matched_markers[label]

    program_scores: dict[str, np.ndarray] = {}
    for program, markers in PROGRAM_MARKERS.items():
        score, matched = panel_score(matrix, gene_index, markers)
        program_scores[program] = score
        summary[f"matched_markers_{program}"] = matched

    masks = {
        "cancer": assignments == "cancer",
        "cytotoxic_t": assignments == "cytotoxic_t",
        "monocyte": assignments == "monocyte",
        "macrophage": assignments == "macrophage",
    }
    program_to_compartment = {
        "tumor_growth": "cancer",
        "tumor_stress": "cancer",
        "tcell_cytotoxicity": "cytotoxic_t",
        "tcell_exhaustion": "cytotoxic_t",
        "monocyte_recruitment": "monocyte",
        "monocyte_antigen_presentation": "monocyte",
        "macrophage_polarization": "macrophage",
        "macrophage_inflammation": "macrophage",
    }
    for program, compartment in program_to_compartment.items():
        score = _mean_over_mask(program_scores[program], masks[compartment])
        if score == 0.0:
            score = float(program_scores[program].mean())
        summary[f"program_{program}"] = score

    summary["interaction_mhci_cross_dressing"] = (
        float(summary["frac_monocyte"])
        * float(summary["frac_cytotoxic_t"])
        * float(summary["program_monocyte_antigen_presentation"])
        * float(summary["program_tcell_cytotoxicity"])
    )
    return summary


def summarize_download_root(
    downloads_root: Path,
    output_csv: Path,
    sample_limit: int | None = None,
) -> list[dict[str, float | str]]:
    sample_dirs = discover_10x_sample_dirs(downloads_root)
    if sample_limit:
        sample_dirs = sample_dirs[:sample_limit]
    if not sample_dirs:
        raise RuntimeError(f"No 10x sample directories found under {downloads_root}")
    summaries = [summarize_10x_sample(sample_dir, downloads_root) for sample_dir in sample_dirs]
    ensure_dir(output_csv.parent)
    write_csv(output_csv, summaries)
    return summaries

