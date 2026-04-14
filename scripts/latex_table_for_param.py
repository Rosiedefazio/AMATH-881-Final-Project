# ------------------------------------------------------------
#  Imports
# ------------------------------------------------------------
from __future__ import annotations

from typing import Any, Dict, Sequence

import pandas as pd


# ------------------------------------------------------------
#  Helper utilities
# ------------------------------------------------------------
def _strip_trailing_a(key: str) -> str:
    """Remove a single trailing 'a' from *key* if present."""
    return key[:-1] if key.endswith("a") else key


def _is_scalar(v: Any) -> bool:
    """True if *v* is a scalar (not an iterable that pandas can use as a column)."""
    return not isinstance(v, (list, tuple, dict, set, pd.Series, pd.Index))


# Mapping from the *cleaned* key to the LaTeX command for the Greek letter.
# Keys that are not Greek symbols are left unchanged.
_GREEK_LATEX_MAP: dict[str, str] = {
    "alpha": r"\alpha",
    "beta": r"\beta",
    "gamma": r"\gamma",
    "delta": r"\delta",
    "epsilon": r"\epsilon",
    "zeta": r"\zeta",
    "eta": r"\eta",
    "theta": r"\theta",
    "iota": r"\iota",
    "kappa": r"\kappa",
    "lam": r"\lambda",   # we receive “lam” after stripping “a” from “lama”
    "mu": r"\mu",
    "nu": r"\nu",
    "xi": r"\xi",
    "omicron": r"o",      # rarely used in math, keep plain “o”
    "pi": r"\pi",
    "rho": r"\rho",
    "sigma": r"\sigma",
    "tau": r"\tau",
    "upsilon": r"\upsilon",
    "phi": r"\phi",
    "chi": r"\chi",
    "psi": r"\psi",
    "omega": r"\omega",
}


def _to_latex_header(name: str) -> str:
    """
    Convert a cleaned column name to the LaTeX header that will appear
    in the table.  If the name is a known Greek symbol, return its
    LaTeX command; otherwise return the name unchanged.
    """
    return _GREEK_LATEX_MAP.get(name, name)


# ------------------------------------------------------------
#  Main conversion function
# ------------------------------------------------------------
def dict_to_latex_greek(
    raw_dict: Dict[str, Any],
    *,
    index: bool = False,
    caption: str | None = None,
    label: str | None = None,
    float_format: str = "%.6f",
) -> str:
    """
    Convert a dictionary to a LaTeX table, automatically
    - stripping a trailing ``a`` from each key,
    - mapping the cleaned key to a LaTeX Greek‑letter command when possible.

    The function works with pure scalars (producing a one‑row table)
    or with iterables (all columns must share the same length).

    Parameters
    ----------
    raw_dict
        Mapping where keys become column names.
    index
        Include the DataFrame index in the LaTeX output?
    caption
        Optional table caption.
    label
        Optional LaTeX label (e.g. ``"tab:my_table"``).
    float_format
        Formatting string for floating‑point numbers.

    Returns
    -------
    str
        LaTeX source for the table.
    """
    if not raw_dict:
        raise ValueError("Input dictionary is empty.")

    # --------------------------------------------------------
    # 1️⃣ Normalise data – scalar → one‑element list, iterable → list
    # --------------------------------------------------------
    if all(_is_scalar(v) for v in raw_dict.values()):
        # All scalars → single‑row table
        data = {k: [v] for k, v in raw_dict.items()}
    else:
        # At least one iterable – ensure all iterables have the same length
        lengths = {len(v) for v in raw_dict.values() if not _is_scalar(v)}
        if len(lengths) != 1:
            raise ValueError(
                "All iterable columns must contain the same number of rows. "
                f"Found lengths: {sorted(lengths)}"
            )
        target_len = next(iter(lengths))
        data = {}
        for k, v in raw_dict.items():
            if _is_scalar(v):
                data[k] = [v] * target_len
            else:
                data[k] = list(v)

    # --------------------------------------------------------
    # 2️⃣ Build the DataFrame
    # --------------------------------------------------------
    df = pd.DataFrame(data)

    # --------------------------------------------------------
    # 3️⃣ Rename columns:
    #     a) strip trailing “a”
    #     b) replace the cleaned name with the LaTeX Greek command when possible
    # --------------------------------------------------------
    rename_map = {
        old: _to_latex_header(_strip_trailing_a(old)) for old in df.columns
    }
    df = df.rename(columns=rename_map)

    # --------------------------------------------------------
    # 4️⃣ Export to LaTeX
    # --------------------------------------------------------
    latex = df.to_latex(
        index=index,
        caption=caption,
        label=label,
        float_format=float_format,
        column_format="l" * len(df.columns),  # left‑align all columns
        escape=False,                      # we already have LaTeX commands
    )
    return latex


# ------------------------------------------------------------
#  Example usage with the dictionary you posted
# ------------------------------------------------------------
if __name__ == "__main__":
    pars = {
        "r": 0.08,
        "b": 0.01,
        "gammaa": 0,
        "alphaa": 0,
        "lama": 0.2,
        "betaa": 0.015900262162975,
        "deltaa": 1,
        "rhoa": 0.1,
        "etaa": 0.1,
        "omegaa": 1,
    }

    latex_table = dict_to_latex_greek(
        pars,
        caption="Model parameters",
        label="tab:model-parameters",
        index=False,
        float_format="%.6f",
    )

    print(latex_table)