from __future__ import annotations

import os
from pathlib import Path

from tumor_atlas_ode.cli import main


def test_run_pipeline_without_downloads_returns_helpful_error(tmp_path: Path, monkeypatch, capsys) -> None:
    publication_id = "demo-publication"
    downloads_dir = tmp_path / "data" / "raw" / publication_id / "downloads"
    metadata_dir = tmp_path / "data" / "raw" / publication_id / "metadata"
    downloads_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    monkeypatch.setattr(
        "tumor_atlas_ode.cli._fetch_and_materialize",
        lambda _publication_id, _root: (
            {
                "downloads_dir": downloads_dir,
                "metadata_dir": metadata_dir,
            },
            [],
        ),
    )

    seen_filters: dict[str, list[str] | None] = {}

    def fake_filter_open_records(records, assay_names=None, levels=None, timepoints=None):
        seen_filters["assay_names"] = assay_names
        seen_filters["levels"] = levels
        seen_filters["timepoints"] = timepoints
        return [
            {
                "assayName": "scRNA-seq",
                "level": "Level 3",
                "TimepointLabel": "Initial_Diagnosis",
                "ParticipantID": "HTA4_1001",
                "Filename": "single_cell_RNAseq_level_3_NBL/SAMPLE_A/matrix.mtx.gz",
                "synapseId": "syn1",
                "downloadSource": "Synapse",
            }
        ]

    monkeypatch.setattr("tumor_atlas_ode.htan.filter_open_records", fake_filter_open_records)

    exit_code = main(["run-pipeline", publication_id, "--root", str(tmp_path)])

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "run-pipeline demo-publication --download" in captured.err
    assert seen_filters["assay_names"] == ["scRNA-seq"]
    assert seen_filters["levels"] == ["Level 3"]
    assert seen_filters["timepoints"] is None


def test_run_pipeline_download_loads_synapse_token_from_dotenv(tmp_path: Path, monkeypatch) -> None:
    publication_id = "demo-publication"
    downloads_dir = tmp_path / "data" / "raw" / publication_id / "downloads"
    metadata_dir = tmp_path / "data" / "raw" / publication_id / "metadata"
    downloads_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)
    (tmp_path / ".env").write_text("SYNAPSE_AUTH_TOKEN=token-from-dotenv\n", encoding="utf-8")

    monkeypatch.delenv("SYNAPSE_AUTH_TOKEN", raising=False)

    monkeypatch.setattr(
        "tumor_atlas_ode.cli._fetch_and_materialize",
        lambda _publication_id, _root: (
            {
                "downloads_dir": downloads_dir,
                "metadata_dir": metadata_dir,
            },
            [],
        ),
    )
    monkeypatch.setattr(
        "tumor_atlas_ode.htan.filter_open_records",
        lambda records, assay_names=None, levels=None, timepoints=None: [
            {
                "assayName": "scRNA-seq",
                "level": "Level 3",
                "TimepointLabel": "Initial_Diagnosis",
                "ParticipantID": "HTA4_1001",
                "Filename": "single_cell_RNAseq_level_3_NBL/SAMPLE_A/matrix.mtx.gz",
                "synapseId": "syn1",
                "downloadSource": "Synapse",
            }
        ],
    )

    seen: dict[str, object] = {}

    def fake_download(records, paths, email, password, auth_token, limit):
        seen["auth_token_arg"] = auth_token
        seen["env_token"] = os.getenv("SYNAPSE_AUTH_TOKEN")
        sample_dir = downloads_dir / "scrna-seq" / "level-3" / "initial-diagnosis" / "participant-1" / "sample-1"
        sample_dir.mkdir(parents=True, exist_ok=True)
        for name in ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]:
            (sample_dir / name).write_text("placeholder", encoding="utf-8")
        return {"downloaded": 3, "failed": 0, "attempted": 3}

    monkeypatch.setattr("tumor_atlas_ode.cli._download_if_requested", fake_download)
    monkeypatch.setattr("tumor_atlas_ode.rna.summarize_download_root", lambda *args, **kwargs: [])
    monkeypatch.setattr("tumor_atlas_ode.ode_priors.load_sample_summaries", lambda path: [])
    monkeypatch.setattr(
        "tumor_atlas_ode.ode_priors.build_ode_priors",
        lambda summaries, group_by="timepoint_label": {"baseline_group": "median-of-groups", "group_by": group_by, "groups": []},
    )
    monkeypatch.setattr("tumor_atlas_ode.ode_priors.write_ode_prior_artifacts", lambda priors, json_path, csv_path: None)

    exit_code = main(["run-pipeline", publication_id, "--root", str(tmp_path), "--download", "--limit", "3"])

    assert exit_code == 0
    assert seen["auth_token_arg"] is None
    assert seen["env_token"] == "token-from-dotenv"
