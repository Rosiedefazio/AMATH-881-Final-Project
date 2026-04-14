from __future__ import annotations

import json
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from tumor_atlas_ode.utils import ensure_dir, slugify, strip_all_suffixes, write_csv, write_json

HTAN_BASE_URL = "https://humantumoratlas.org"
USER_AGENT = "tumor-atlas-ode/0.1 (+https://humantumoratlas.org)"


@dataclass
class PublicationBundle:
    publication_id: str
    build_id: str
    page_props: dict[str, Any]
    payload: dict[str, Any]


def fetch_text(url: str, timeout: int = 60) -> str:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(request, timeout=timeout) as response:
            return response.read().decode("utf-8")
    except HTTPError as exc:
        raise RuntimeError(f"HTTP error fetching {url}: {exc.code}") from exc
    except URLError as exc:
        raise RuntimeError(f"Network error fetching {url}: {exc}") from exc


def extract_build_id(html: str) -> str:
    match = re.search(r"/_next/static/([^/]+)/_buildManifest\.js", html)
    if not match:
        raise RuntimeError("Could not extract HTAN Next.js build id from publication page.")
    return match.group(1)


def publication_page_url(publication_id: str) -> str:
    return f"{HTAN_BASE_URL}/publications/{publication_id}?tab=overview"


def publication_json_url(publication_id: str, build_id: str) -> str:
    return f"{HTAN_BASE_URL}/_next/data/{build_id}/publications/{publication_id}.json?tab=overview"


def fetch_publication_bundle(publication_id: str) -> PublicationBundle:
    html = fetch_text(publication_page_url(publication_id))
    build_id = extract_build_id(html)
    payload = json.loads(fetch_text(publication_json_url(publication_id, build_id)))
    page_props = payload["pageProps"]
    return PublicationBundle(
        publication_id=publication_id,
        build_id=build_id,
        page_props=page_props,
        payload=payload,
    )


def build_specimen_index(page_props: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {row["BiospecimenID"]: row for row in page_props.get("specimen", []) if row.get("BiospecimenID")}


def infer_participant_id(record: dict[str, Any]) -> str:
    demographics_ids = record.get("demographicsIds") or []
    if demographics_ids:
        return demographics_ids[0]
    if record.get("ParticipantID"):
        return str(record["ParticipantID"])
    return "unknown-participant"


def infer_timepoint_label(record: dict[str, Any], specimen_index: dict[str, dict[str, Any]]) -> str:
    labels = []
    for biospecimen_id in record.get("biospecimenIds") or []:
        specimen = specimen_index.get(biospecimen_id, {})
        label = specimen.get("TimepointLabel")
        if label and label not in labels:
            labels.append(label)
    if not labels:
        return "Unknown"
    if len(labels) == 1:
        return labels[0]
    return "+".join(labels)


def infer_sample_key(record: dict[str, Any]) -> str:
    filename = record.get("Filename") or ""
    path = PurePosixPath(filename)
    generic_parents = {
        "single_cell_RNAseq_level_3_NBL",
        "single_cell_RNAseq_level_4_NBL",
        "seurat_objects_count_based",
        "seurat_objects_sct_based",
        "rds",
    }
    if path.name in {"features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz"} and len(path.parts) >= 2:
        return path.parent.name
    if len(path.parts) >= 2 and path.parent.name not in generic_parents:
        return path.parent.name
    return strip_all_suffixes(path.name or record.get("DataFileID", "sample"))


def infer_download_subdir(record: dict[str, Any]) -> Path:
    return (
        Path(slugify(record.get("assayName", "unknown-assay")))
        / slugify(record.get("level", "unknown-level"))
        / slugify(record.get("TimepointLabel", "unknown-timepoint"))
        / slugify(record.get("ParticipantID", "unknown-participant"))
        / slugify(record.get("SampleKey", "unknown-sample"))
    )


def annotate_assay_records(page_props: dict[str, Any]) -> list[dict[str, Any]]:
    specimen_index = build_specimen_index(page_props)
    annotated: list[dict[str, Any]] = []
    for record in page_props.get("assays", []):
        item = dict(record)
        item["ParticipantID"] = infer_participant_id(item)
        item["TimepointLabel"] = infer_timepoint_label(item, specimen_index)
        item["SampleKey"] = infer_sample_key(item)
        item["DownloadSubdir"] = str(infer_download_subdir(item))
        annotated.append(item)
    return annotated


def filter_open_records(
    records: list[dict[str, Any]],
    assay_names: list[str] | None = None,
    levels: list[str] | None = None,
    timepoints: list[str] | None = None,
) -> list[dict[str, Any]]:
    filtered = [row for row in records if row.get("downloadSource") == "Synapse" and row.get("synapseId")]
    if assay_names:
        allowed = set(assay_names)
        filtered = [row for row in filtered if row.get("assayName") in allowed]
    if levels:
        allowed = set(levels)
        filtered = [row for row in filtered if row.get("level") in allowed]
    if timepoints:
        allowed = set(timepoints)
        filtered = [row for row in filtered if row.get("TimepointLabel") in allowed]
    return sorted(filtered, key=lambda row: (
        row.get("assayName", ""),
        row.get("level", ""),
        row.get("TimepointLabel", ""),
        row.get("ParticipantID", ""),
        row.get("Filename", ""),
    ))


def group_records_for_manifests(records: list[dict[str, Any]]) -> dict[tuple[str, str, str], list[dict[str, Any]]]:
    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        key = (
            slugify(record.get("assayName", "unknown-assay")),
            slugify(record.get("level", "unknown-level")),
            slugify(record.get("TimepointLabel", "unknown-timepoint")),
        )
        grouped[key].append(record)
    return grouped


def materialize_publication(root: Path, bundle: PublicationBundle) -> dict[str, Path]:
    publication_root = ensure_dir(root / "data" / "raw" / bundle.publication_id)
    metadata_dir = ensure_dir(publication_root / "metadata")
    manifests_dir = ensure_dir(publication_root / "manifests")
    downloads_dir = ensure_dir(publication_root / "downloads")

    annotated_assays = annotate_assay_records(bundle.page_props)
    open_assays = filter_open_records(annotated_assays)

    write_json(metadata_dir / "publication.json", bundle.payload)
    write_json(
        metadata_dir / "fetch_metadata.json",
        {
            "publication_id": bundle.publication_id,
            "build_id": bundle.build_id,
            "publication_json_url": publication_json_url(bundle.publication_id, bundle.build_id),
        },
    )
    write_csv(metadata_dir / "assays.csv", annotated_assays)
    write_csv(metadata_dir / "open_access_assays.csv", open_assays)
    write_csv(metadata_dir / "specimen.csv", bundle.page_props.get("specimen", []))
    write_csv(metadata_dir / "cases.csv", bundle.page_props.get("cases", []))

    for (assay_slug, level_slug, timepoint_slug), rows in group_records_for_manifests(annotated_assays).items():
        write_csv(manifests_dir / assay_slug / level_slug / timepoint_slug / "manifest.csv", rows)

    return {
        "publication_root": publication_root,
        "metadata_dir": metadata_dir,
        "manifests_dir": manifests_dir,
        "downloads_dir": downloads_dir,
    }
