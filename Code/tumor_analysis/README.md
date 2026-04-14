# Tumor Atlas ODE

Python repo for two linked tasks:

1. Pull HTAN publication metadata and open-access scRNA-seq assets from the CHOP neuroblastoma publication page.
2. Convert downloaded scRNA-seq timepoint data into relative initial conditions and parameter scales for a mechanistic ODE model with:
   - cancer volume
   - cytotoxic T-cell volume
   - monocyte volume
   - macrophage volume
   - cross-interactions, including a monocyte/CD8 interaction proxy for MHC-I cross-dressing

## What It Uses

- HTAN publication metadata is fetched from the publication's Next.js JSON endpoint.
- Open-access file downloads use `synapseId` values exposed by HTAN records.
- Level 3 scRNA-seq 10x matrices are summarized with marker-based compartment scoring.
- ODE priors are generated as relative scales, anchored to `Initial_Diagnosis` when that timepoint exists.

## Run Locally With venv

If you have Python 3.12+ available locally, you can run the scraper and analysis without Docker:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[synapse,dev]'
tumor-atlas-ode --help
```

If `python3.12` on your shell still points to the system Python, use the full path to your configured interpreter instead, for example:

```bash
~/.pyenv/shims/python3.12 -m venv .venv
```

Once the environment is active, every CLI example below can be run directly by replacing the Docker wrapper with `tumor-atlas-ode`.

## Run With Docker

Build the image:

```bash
docker build -t tumor-atlas-ode .
```

Run the CLI against the current repo checkout:

```bash
docker run --rm -it \
  -v "$PWD:/app" \
  -w /app \
  tumor-atlas-ode --help
```

The image installs the package with the `synapse` extra, so the download commands are available inside the container.

## Synapse Credentials

For local venv usage, export:

```bash
export SYNAPSE_AUTH_TOKEN='...'
```

Local CLI runs also auto-load `./.env` from the selected `--root` directory if that file exists, without overriding already-exported shell variables.

For Docker usage, pass:

```bash
--env SYNAPSE_AUTH_TOKEN='...'
```

You can also use `--env-file .env`, or rely on a configured local Synapse profile.

Note: metadata may be visible anonymously while file bytes still require authenticated download permission.

## CLI

Fetch HTAN metadata and build manifests:

```bash
docker run --rm -it \
  -v "$PWD:/app" \
  -w /app \
  tumor-atlas-ode \
  fetch-publication hta4_2025_nature-genetics_wenbao-yu
```

Download open-access scRNA-seq Level 3 files from Synapse:

```bash
docker run --rm -it \
  -v "$PWD:/app" \
  -w /app \
  --env-file .env \
  tumor-atlas-ode \
  download-open hta4_2025_nature-genetics_wenbao-yu \
  --assay scRNA-seq \
  --level "Level 3"
```

Summarize downloaded Level 3 matrices:

```bash
docker run --rm -it \
  -v "$PWD:/app" \
  -w /app \
  tumor-atlas-ode \
  summarize-scrna --publication-id hta4_2025_nature-genetics_wenbao-yu
```

Build ODE priors from sample summaries:

```bash
docker run --rm -it \
  -v "$PWD:/app" \
  -w /app \
  tumor-atlas-ode \
  build-ode-priors \
  --summary-csv data/processed/hta4_2025_nature-genetics_wenbao-yu/rna/sample_summaries.csv
```

End-to-end convenience command:

```bash
docker run --rm -it \
  -v "$PWD:/app" \
  -w /app \
  --env-file .env \
  tumor-atlas-ode \
  run-pipeline hta4_2025_nature-genetics_wenbao-yu \
  --download \
  --level "Level 3"
```

## Output Layout

```text
data/
  raw/
    <publication_id>/
      metadata/
        publication.json
        assays.csv
        cases.csv
        specimen.csv
        open_access_assays.csv
        downloaded_files.csv
      manifests/
        <assay>/<level>/<timepoint>/manifest.csv
      downloads/
        <assay>/<level>/<timepoint>/<participant>/<sample_key>/*
  processed/
    <publication_id>/
      rna/
        sample_summaries.csv
      ode/
        priors_by_timepoint.json
        priors_by_timepoint.csv
```

## Modeling Assumptions

- Initial conditions are relative compartment volumes, not absolute mm^3 values.
- Parameter outputs are relative scales intended to initialize or constrain a downstream ODE calibration step.
- MHC-I cross-dressing is represented by a proxy interaction score:
  - monocyte abundance
  - monocyte antigen-presentation program
  - CD8 abundance
  - CD8 cytotoxic activation program

## Current Limitation

The scRNA summarizer is marker-based and intentionally lightweight. If you already have curated cell-type annotations or Seurat/Scanpy outputs, the prior builder can be extended to consume those directly and replace the heuristic classifier.

## Local Verification

```bash
source .venv/bin/activate
python -m pytest -q
tumor-atlas-ode fetch-publication hta4_2025_nature-genetics_wenbao-yu --root .
```
