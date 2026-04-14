from __future__ import annotations

import csv
import json
import os
import re
from pathlib import Path
from typing import Any, Iterable


def slugify(value: str) -> str:
    text = value.strip().lower()
    text = re.sub(r"[^a-z0-9]+", "-", text)
    return text.strip("-") or "unknown"


def strip_all_suffixes(name: str) -> str:
    path = Path(name)
    stem = path.name
    while True:
        next_stem = Path(stem).stem
        if next_stem == stem:
            return stem
        stem = next_stem


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def json_string(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list)):
        return json.dumps(value, sort_keys=True)
    return str(value)


def flatten_record(record: dict[str, Any], prefix: str = "") -> dict[str, str]:
    flat: dict[str, str] = {}
    for key, value in record.items():
        next_key = f"{prefix}__{key}" if prefix else key
        if isinstance(value, dict):
            flat.update(flatten_record(value, prefix=next_key))
        else:
            flat[next_key] = json_string(value)
    return flat


def write_json(path: Path, payload: Any) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    row_list = [flatten_record(row) for row in rows]
    ensure_dir(path.parent)
    if not row_list:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = sorted({key for row in row_list for key in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(row_list)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader)


def parse_jsonish(value: str) -> Any:
    if value is None or value == "":
        return value
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        return value


def try_float(value: Any) -> float:
    if value in ("", None):
        return 0.0
    if isinstance(value, (int, float)):
        return float(value)
    return float(str(value))


def load_env_file(path: Path, override: bool = False) -> dict[str, str]:
    loaded: dict[str, str] = {}
    if not path.exists():
        return loaded

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[7:].strip()
        if "=" not in line:
            continue

        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key:
            continue

        if value and value[0] == value[-1] and value[0] in {"'", '"'}:
            value = value[1:-1]
        if not override and key in os.environ:
            continue

        os.environ[key] = value
        loaded[key] = value

    return loaded
