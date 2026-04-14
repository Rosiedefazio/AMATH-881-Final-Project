from __future__ import annotations

import os
import shutil
from pathlib import Path

from tumor_atlas_ode.htan import infer_download_subdir
from tumor_atlas_ode.utils import ensure_dir, write_csv


def _load_synapse_client(email: str | None, password: str | None, auth_token: str | None):
    try:
        import synapseclient
    except ImportError as exc:
        raise RuntimeError(
            "synapseclient is not installed. Run `pip install -e '.[synapse]'` to enable downloads."
        ) from exc

    syn = synapseclient.Synapse()
    if not any([email, password, auth_token]):
        return syn

    login_kwargs = {
        "silent": True,
    }
    if auth_token:
        login_kwargs["authToken"] = auth_token
    if email:
        login_kwargs["email"] = email
    if password and not auth_token:
        raise RuntimeError(
            "This synapseclient version does not support password login in the API wrapper here. "
            "Use SYNAPSE_AUTH_TOKEN or a configured Synapse profile instead."
        )
    syn.login(**login_kwargs)
    return syn


def download_synapse_records(
    records: list[dict[str, str]],
    downloads_root: Path,
    metadata_dir: Path,
    email: str | None = None,
    password: str | None = None,
    auth_token: str | None = None,
    limit: int | None = None,
) -> dict[str, int]:
    email = email or os.getenv("SYNAPSE_EMAIL")
    password = password or os.getenv("SYNAPSE_PASSWORD")
    auth_token = auth_token or os.getenv("SYNAPSE_AUTH_TOKEN")
    syn = _load_synapse_client(email=email, password=password, auth_token=auth_token)

    results: list[dict[str, str]] = []
    selected = records[:limit] if limit else records
    success_count = 0
    failure_count = 0

    for record in selected:
        synapse_id = record.get("synapseId", "")
        target_dir = ensure_dir(downloads_root / infer_download_subdir(record))
        expected_name = Path(record.get("Filename", "")).name or f"{record.get('DataFileID', synapse_id)}"
        try:
            entity = syn.get(synapse_id, downloadLocation=str(target_dir), ifcollision="overwrite.local")
            if not getattr(entity, "path", None):
                raise RuntimeError(
                    f"Synapse entity {synapse_id} is visible but was not downloaded. "
                    "This usually means the file requires authenticated download permission."
                )
            downloaded_path = Path(entity.path)
            expected_path = target_dir / expected_name
            if downloaded_path.name != expected_name:
                if expected_path.exists():
                    expected_path.unlink()
                shutil.move(str(downloaded_path), str(expected_path))
                downloaded_path = expected_path
            results.append(
                {
                    **record,
                    "local_path": str(downloaded_path),
                    "download_status": "downloaded",
                }
            )
            success_count += 1
        except Exception as exc:  # noqa: BLE001
            results.append(
                {
                    **record,
                    "local_path": "",
                    "download_status": f"error: {exc}",
                }
            )
            failure_count += 1

    write_csv(metadata_dir / "downloaded_files.csv", results)
    return {"downloaded": success_count, "failed": failure_count, "attempted": len(selected)}
