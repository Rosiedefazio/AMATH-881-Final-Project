from __future__ import annotations

import argparse
import sys
from pathlib import Path

DEFAULT_PIPELINE_ASSAYS = ["scRNA-seq"]
DEFAULT_PIPELINE_LEVELS = ["Level 3"]


def add_common_filters(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--assay", action="append", help="Assay name filter, repeatable.")
    parser.add_argument("--level", action="append", help="Level filter, repeatable.")
    parser.add_argument("--timepoint", action="append", help="Timepoint filter, repeatable.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="HTAN ingestion and ODE prior generation.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    fetch_parser = subparsers.add_parser("fetch-publication", help="Fetch HTAN publication metadata.")
    fetch_parser.add_argument("publication_id")
    fetch_parser.add_argument("--root", default=".", help="Repo root for data output.")

    download_parser = subparsers.add_parser("download-open", help="Download Synapse-backed open-access files.")
    download_parser.add_argument("publication_id")
    download_parser.add_argument("--root", default=".", help="Repo root for data output.")
    download_parser.add_argument("--synapse-email", default=None)
    download_parser.add_argument("--synapse-password", default=None)
    download_parser.add_argument("--synapse-token", default=None)
    download_parser.add_argument("--limit", type=int, default=None)
    add_common_filters(download_parser)

    summarize_parser = subparsers.add_parser("summarize-scrna", help="Summarize downloaded Level 3 10x data.")
    summarize_parser.add_argument("--publication-id", required=True)
    summarize_parser.add_argument("--root", default=".", help="Repo root for data output.")
    summarize_parser.add_argument("--sample-limit", type=int, default=None)

    priors_parser = subparsers.add_parser("build-ode-priors", help="Build ODE priors from sample summaries.")
    priors_parser.add_argument("--summary-csv", required=True)
    priors_parser.add_argument("--output-json", default=None)
    priors_parser.add_argument("--output-csv", default=None)
    priors_parser.add_argument("--group-by", default="timepoint_label")

    pipeline_parser = subparsers.add_parser(
        "run-pipeline",
        help="Fetch metadata, optionally download supported scRNA-seq inputs, summarize, and build priors.",
    )
    pipeline_parser.add_argument("publication_id")
    pipeline_parser.add_argument("--root", default=".", help="Repo root for data output.")
    pipeline_parser.add_argument(
        "--download",
        action="store_true",
        help="Download open scRNA-seq files before summarizing. Fresh runs need this unless downloads already exist.",
    )
    pipeline_parser.add_argument("--synapse-email", default=None)
    pipeline_parser.add_argument("--synapse-password", default=None)
    pipeline_parser.add_argument("--synapse-token", default=None)
    pipeline_parser.add_argument("--limit", type=int, default=None)
    pipeline_parser.add_argument("--sample-limit", type=int, default=None)
    add_common_filters(pipeline_parser)

    return parser


def _fetch_and_materialize(publication_id: str, root: Path) -> tuple[dict[str, Path], list[dict[str, str]]]:
    from tumor_atlas_ode.htan import (
        annotate_assay_records,
        fetch_publication_bundle,
        materialize_publication,
    )

    bundle = fetch_publication_bundle(publication_id)
    paths = materialize_publication(root, bundle)
    assays = annotate_assay_records(bundle.page_props)
    return paths, assays


def _load_root_env(root: Path) -> None:
    from tumor_atlas_ode.utils import load_env_file

    load_env_file(root / ".env", override=False)


def _download_if_requested(
    records: list[dict[str, str]],
    paths: dict[str, Path],
    email: str | None,
    password: str | None,
    auth_token: str | None,
    limit: int | None,
) -> dict[str, int]:
    from tumor_atlas_ode.synapse_download import download_synapse_records

    return download_synapse_records(
        records=records,
        downloads_root=paths["downloads_dir"],
        metadata_dir=paths["metadata_dir"],
        email=email,
        password=password,
        auth_token=auth_token,
        limit=limit,
    )


def _pipeline_filters(args: argparse.Namespace) -> tuple[list[str], list[str], list[str] | None]:
    assay_names = args.assay or list(DEFAULT_PIPELINE_ASSAYS)
    levels = args.level or list(DEFAULT_PIPELINE_LEVELS)
    return assay_names, levels, args.timepoint


def _ensure_pipeline_inputs(
    downloads_root: Path,
    publication_id: str,
    attempted_download: bool,
    limit: int | None,
) -> None:
    from tumor_atlas_ode.rna import discover_10x_sample_dirs

    if discover_10x_sample_dirs(downloads_root):
        return
    if attempted_download:
        message = (
            f"No complete 10x sample directories were found under {downloads_root}. "
            "The RNA summarizer requires matching `barcodes.tsv.gz`, `features.tsv.gz`, "
            "and `matrix.mtx.gz` files in the same sample directory. "
        )
        if limit is not None and limit % 3 != 0:
            message += (
                "You used `--limit`; for Level 3 10x inputs, use a multiple of 3 such as `--limit 3` "
                "to fetch one complete sample. "
            )
        message += "Check your download filters and Synapse permissions."
        raise RuntimeError(message)
    raise RuntimeError(
        f"No downloaded Level 3 scRNA-seq 10x sample directories were found under {downloads_root}. "
        f"Run `tumor-atlas-ode run-pipeline {publication_id} --download` to fetch them first, "
        "or pre-download them with `tumor-atlas-ode download-open ... --assay scRNA-seq --level \"Level 3\"`."
    )


def command_fetch_publication(args: argparse.Namespace) -> int:
    root = Path(args.root).resolve()
    paths, assays = _fetch_and_materialize(args.publication_id, root)
    print(f"Fetched metadata for {args.publication_id}")
    print(f"Metadata directory: {paths['metadata_dir']}")
    print(f"Assay rows: {len(assays)}")
    return 0


def command_download_open(args: argparse.Namespace) -> int:
    from tumor_atlas_ode.htan import filter_open_records

    root = Path(args.root).resolve()
    _load_root_env(root)
    paths, assays = _fetch_and_materialize(args.publication_id, root)
    records = filter_open_records(assays, assay_names=args.assay, levels=args.level, timepoints=args.timepoint)
    stats = _download_if_requested(
        records=records,
        paths=paths,
        email=args.synapse_email,
        password=args.synapse_password,
        auth_token=args.synapse_token,
        limit=args.limit,
    )
    print(f"Downloaded {stats['downloaded']} files, failed {stats['failed']}, attempted {stats['attempted']}")
    print(f"Download log: {paths['metadata_dir'] / 'downloaded_files.csv'}")
    return 1 if stats["failed"] else 0


def command_summarize_scrna(args: argparse.Namespace) -> int:
    from tumor_atlas_ode.rna import summarize_download_root

    root = Path(args.root).resolve()
    publication_root = root / "data" / "raw" / args.publication_id
    downloads_root = publication_root / "downloads"
    output_csv = root / "data" / "processed" / args.publication_id / "rna" / "sample_summaries.csv"
    summaries = summarize_download_root(downloads_root=downloads_root, output_csv=output_csv, sample_limit=args.sample_limit)
    print(f"Wrote {len(summaries)} sample summaries to {output_csv}")
    return 0


def command_build_ode_priors(args: argparse.Namespace) -> int:
    from tumor_atlas_ode.ode_priors import build_ode_priors, load_sample_summaries, write_ode_prior_artifacts

    summary_csv = Path(args.summary_csv).resolve()
    publication_id = summary_csv.parent.parent.name
    output_json = Path(args.output_json).resolve() if args.output_json else summary_csv.parent.parent / "ode" / "priors_by_timepoint.json"
    output_csv = Path(args.output_csv).resolve() if args.output_csv else summary_csv.parent.parent / "ode" / "priors_by_timepoint.csv"
    priors = build_ode_priors(load_sample_summaries(summary_csv), group_by=args.group_by)
    write_ode_prior_artifacts(priors, json_path=output_json, csv_path=output_csv)
    print(f"Wrote ODE priors for {publication_id}")
    print(f"JSON: {output_json}")
    print(f"CSV: {output_csv}")
    return 0


def command_run_pipeline(args: argparse.Namespace) -> int:
    from tumor_atlas_ode.htan import filter_open_records
    from tumor_atlas_ode.ode_priors import build_ode_priors, load_sample_summaries, write_ode_prior_artifacts
    from tumor_atlas_ode.rna import summarize_download_root

    root = Path(args.root).resolve()
    _load_root_env(root)
    paths, assays = _fetch_and_materialize(args.publication_id, root)
    assay_names, levels, timepoints = _pipeline_filters(args)
    records = filter_open_records(assays, assay_names=assay_names, levels=levels, timepoints=timepoints)
    if args.download:
        if not records:
            raise RuntimeError(
                "No open Synapse-backed scRNA-seq Level 3 records matched the pipeline filters. "
                "Adjust `--assay`, `--level`, or `--timepoint` and try again."
            )
        stats = _download_if_requested(
            records=records,
            paths=paths,
            email=args.synapse_email,
            password=args.synapse_password,
            auth_token=args.synapse_token,
            limit=args.limit,
        )
        print(f"Downloaded {stats['downloaded']} files, failed {stats['failed']}, attempted {stats['attempted']}")
        if stats["failed"]:
            return 1
    _ensure_pipeline_inputs(
        downloads_root=paths["downloads_dir"],
        publication_id=args.publication_id,
        attempted_download=args.download,
        limit=args.limit,
    )
    output_csv = root / "data" / "processed" / args.publication_id / "rna" / "sample_summaries.csv"
    summarize_download_root(paths["downloads_dir"], output_csv=output_csv, sample_limit=args.sample_limit)
    priors = build_ode_priors(load_sample_summaries(output_csv), group_by="timepoint_label")
    write_ode_prior_artifacts(
        priors,
        json_path=root / "data" / "processed" / args.publication_id / "ode" / "priors_by_timepoint.json",
        csv_path=root / "data" / "processed" / args.publication_id / "ode" / "priors_by_timepoint.csv",
    )
    print(f"Pipeline complete for {args.publication_id}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "fetch-publication":
            return command_fetch_publication(args)
        if args.command == "download-open":
            return command_download_open(args)
        if args.command == "summarize-scrna":
            return command_summarize_scrna(args)
        if args.command == "build-ode-priors":
            return command_build_ode_priors(args)
        if args.command == "run-pipeline":
            return command_run_pipeline(args)
        parser.error(f"Unknown command: {args.command}")
        return 2
    except RuntimeError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
