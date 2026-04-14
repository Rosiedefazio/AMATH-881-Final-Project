from __future__ import annotations

import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

from tumor_atlas_ode.utils import read_csv, try_float, write_csv, write_json


def _mean(rows: list[dict[str, Any]], field: str) -> float:
    values = [try_float(row.get(field, 0.0)) for row in rows]
    return statistics.fmean(values) if values else 0.0


def _baseline_group_name(group_names: list[str]) -> str | None:
    for name in group_names:
        normalized = name.strip().lower().replace("_", "-").replace(" ", "-")
        if normalized == "initial-diagnosis":
            return name
    return None


def _raw_parameter_bundle(rows: list[dict[str, Any]]) -> dict[str, float]:
    cancer_volume = _mean(rows, "frac_cancer")
    cytotoxic_t_volume = _mean(rows, "frac_cytotoxic_t")
    monocyte_volume = _mean(rows, "frac_monocyte")
    macrophage_volume = _mean(rows, "frac_macrophage")

    tumor_growth = _mean(rows, "program_tumor_growth")
    tumor_stress = _mean(rows, "program_tumor_stress")
    tcell_cytotoxicity = _mean(rows, "program_tcell_cytotoxicity")
    tcell_exhaustion = _mean(rows, "program_tcell_exhaustion")
    monocyte_recruitment = _mean(rows, "program_monocyte_recruitment")
    monocyte_antigen_presentation = _mean(rows, "program_monocyte_antigen_presentation")
    macrophage_polarization = _mean(rows, "program_macrophage_polarization")
    macrophage_inflammation = _mean(rows, "program_macrophage_inflammation")
    mhci_cross_dressing = _mean(rows, "interaction_mhci_cross_dressing")
    eps = 1e-6

    return {
        "cancer_growth_rate": tumor_growth * (1.0 + cancer_volume),
        "cancer_death_rate": tumor_stress / (tumor_growth + eps),
        "cytotoxic_t_growth_rate": tcell_cytotoxicity * (1.0 + cytotoxic_t_volume),
        "cytotoxic_t_death_rate": tcell_exhaustion + eps,
        "monocyte_growth_rate": monocyte_recruitment * (1.0 + monocyte_volume),
        "macrophage_growth_rate": macrophage_polarization * (1.0 + macrophage_volume),
        "macrophage_death_rate": macrophage_inflammation / (macrophage_polarization + eps),
        "tumor_tcell_interaction": tcell_cytotoxicity * cytotoxic_t_volume / (cancer_volume + eps),
        "tumor_monocyte_interaction": monocyte_recruitment * monocyte_volume,
        "tumor_macrophage_interaction": macrophage_polarization * macrophage_volume,
        "monocyte_to_macrophage_transition": monocyte_volume * macrophage_polarization,
        "mhci_cross_dressing_interaction": mhci_cross_dressing
        + (monocyte_antigen_presentation * tcell_cytotoxicity * monocyte_volume * cytotoxic_t_volume),
    }


def build_ode_priors(
    summary_rows: list[dict[str, Any]],
    group_by: str = "timepoint_label",
) -> dict[str, Any]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in summary_rows:
        grouped[str(row.get(group_by, "unknown"))].append(row)

    group_names = sorted(grouped)
    baseline_name = _baseline_group_name(group_names)
    raw_by_group = {group: _raw_parameter_bundle(rows) for group, rows in grouped.items()}

    if baseline_name is not None:
        baseline_metrics = raw_by_group[baseline_name]
    else:
        baseline_metrics = {}
        metric_names = next(iter(raw_by_group.values())).keys()
        for metric_name in metric_names:
            positive = [metrics[metric_name] for metrics in raw_by_group.values() if metrics[metric_name] > 0]
            baseline_metrics[metric_name] = statistics.median(positive) if positive else 1.0

    groups = []
    for group in group_names:
        rows = grouped[group]
        initial_conditions = {
            "cancer_volume": _mean(rows, "frac_cancer"),
            "cytotoxic_t_volume": _mean(rows, "frac_cytotoxic_t"),
            "monocyte_volume": _mean(rows, "frac_monocyte"),
            "macrophage_volume": _mean(rows, "frac_macrophage"),
        }
        raw_parameters = raw_by_group[group]
        parameter_scales = {}
        for metric_name, raw_value in raw_parameters.items():
            baseline = baseline_metrics.get(metric_name, 1.0) or 1.0
            parameter_scales[metric_name] = raw_value / baseline
        groups.append(
            {
                "group": group,
                "sample_count": len(rows),
                "initial_conditions_relative_volume": initial_conditions,
                "parameter_scales": parameter_scales,
                "raw_parameter_scores": raw_parameters,
                "evidence": {
                    "mean_tumor_growth_program": _mean(rows, "program_tumor_growth"),
                    "mean_tcell_cytotoxicity_program": _mean(rows, "program_tcell_cytotoxicity"),
                    "mean_monocyte_antigen_presentation_program": _mean(rows, "program_monocyte_antigen_presentation"),
                    "mean_mhci_cross_dressing_proxy": _mean(rows, "interaction_mhci_cross_dressing"),
                },
            }
        )

    return {
        "group_by": group_by,
        "baseline_group": baseline_name or "median-of-groups",
        "groups": groups,
    }


def load_sample_summaries(path: Path) -> list[dict[str, Any]]:
    return read_csv(path)


def write_ode_prior_artifacts(
    priors: dict[str, Any],
    json_path: Path,
    csv_path: Path,
) -> None:
    write_json(json_path, priors)
    csv_rows = []
    for group in priors["groups"]:
        row = {
            "group": group["group"],
            "sample_count": group["sample_count"],
            **group["initial_conditions_relative_volume"],
            **group["parameter_scales"],
        }
        csv_rows.append(row)
    write_csv(csv_path, csv_rows)
