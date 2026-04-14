"""
Template: Codim-1 and Codim-2 bifurcation analysis for 3D ODEs
using PyDSTool (PyDSCont).

Workflow
--------
  EP-C  (1 free par)   →  finds LP, H, BP on equilibrium branch
  LC-C  (1 free par)   →  limit-cycle branch from first Hopf point
  LP-C  (2 free pars)  →  fold locus; detects Cusp, BT, ZH         [Codim-2]
  H-C   (2 free pars)  →  Hopf locus; detects BT, GH, ZH, HH       [Codim-2]

Edit only sections 1 and 2 for a new model.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter
from typing import Dict, List, Optional, Any

from PyDSTool import *

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1 – SYSTEM DEFINITION  (only part you change per model)
# ─────────────────────────────────────────────────────────────────────────────

SYSTEM = dict(
    name="Lorenz",

    # Parameter name → initial value
    pars={
        "sigma": 10.0,
        "rho":    0.5,   # primary bifurcation parameter (start away from chaos)
        "beta":   2.667,
    },

    # Variable name → RHS string  (PyDSTool math syntax)
    varspecs={
        "x": "sigma * (y - x)",
        "y": "x * (rho - z) - y",
        "z": "x * y - beta * z",
    },

    # Initial condition (used as the seed equilibrium)
    ics={"x": 0.1, "y": 0.1, "z": 0.1},

    # Integration window (used internally by PyDSTool)
    tdata=[0.0, 100.0],
)

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2 – CONTINUATION SETTINGS
# ─────────────────────────────────────────────────────────────────────────────

FREE_PAR_1 = "rho"     # Codim-1 sweep
FREE_PAR_2 = "sigma"   # second free parameter for Codim-2 loci

PAR_DOMAIN = {
    FREE_PAR_1: [0.0, 30.0],
    FREE_PAR_2: [1.0, 20.0],
}

BASE_CONT = dict(
    StepSize    = 1e-2,
    MaxStepSize = 5e-1,
    MaxNumPoints= 300,
    verbosity   = 2,
    SaveEigen   = True,
    LocBifPoints= "all",
    StopAtPoints= ["B"],
)

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3 – HELPERS  (do not edit)
# ─────────────────────────────────────────────────────────────────────────────

def build_ode(sys_dict: dict) -> Any:
    """Construct a PyDSTool Vode_ODEsystem from a plain dict."""
    d = args(name=sys_dict["name"])
    d.pars      = dict(sys_dict["pars"])
    d.varspecs  = dict(sys_dict["varspecs"])
    d.ics       = dict(sys_dict["ics"])
    d.tdata     = list(sys_dict["tdata"])
    # Allow parameters to roam freely; tighten per-model if needed
    d.pdomain   = {k: [-1e9, 1e9] for k in sys_dict["pars"]}
    return Generator.Vode_ODEsystem(d)


def _get_sp(PC, curve: str, label: str) -> Any:
    """Return a special point or None (never raises)."""
    try:
        return PC[curve].getSpecialPoint(label)
    except Exception:
        return None


def run_curve(
    PC,
    name: str,
    curve_type: str,
    freepars: List[str],
    initpoint: Any,
    initdirec: Any = None,
    extra: Optional[Dict] = None,
    overrides: Optional[Dict] = None,
    bidirectional: bool = True,
) -> bool:
    """
    Register, run (forward + optional backward), and clean a continuation curve.

    Returns True on success, False if an exception is raised.
    """
    s = dict(BASE_CONT)
    if overrides:
        s.update(overrides)

    a = args(name=name, type=curve_type)
    a.freepars = freepars
    a.initpoint = initpoint
    for k, v in s.items():
        setattr(a, k, v)
    if initdirec is not None:
        a.initdirec = initdirec
    if extra:
        for k, v in extra.items():
            setattr(a, k, v)

    try:
        PC.newCurve(a)
        t0 = perf_counter()
        PC[name].forward()
        if bidirectional:
            PC[name].backward()
            PC[name].cleanLabels()
        print("  ✓ %s  (%.3f s)" % (name, perf_counter() - t0))
        return True
    except Exception as e:
        print("  ✗ %s skipped: %s" % (name, e))
        return False


def orient(PC, curve: str, sp_label: str, parname: str, increasing: bool = True) -> Optional[Dict]:
    """
    Return the tangent direction dict for a new branch starting at a
    special point, oriented so *parname* increases (or decreases).
    Returns None if the point is unavailable.
    """
    sp = _get_sp(PC, curve, sp_label)
    if sp is None:
        return None
    bif_type = sp_label.rstrip("0123456789")
    try:
        branch = sp.labels[bif_type]["data"].branch
        sign = 1 if increasing else -1
        if branch[parname] * sign < 0:
            return {k: -v for k, v in branch.items()}
        return branch
    except Exception as e:
        print("  [orient] %s:%s – %s" % (curve, sp_label, e))
        return None


def list_special_points(PC, curve: str) -> List[str]:
    """Return sorted list of special-point labels on a curve."""
    try:
        pts = PC[curve].getSpecialPoint()      # returns dict-like object
        return sorted(pts.keys())
    except Exception:
        return []


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4 – BUILD SYSTEM
# ─────────────────────────────────────────────────────────────────────────────

ode = build_ode(SYSTEM)
PC  = ContClass(ode)

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5 – CODIM-1: EQUILIBRIUM BRANCH  (EP-C)
# ─────────────────────────────────────────────────────────────────────────────

print("\n=== EP1: equilibrium continuation ===")
run_curve(
    PC, "EP1", "EP-C",
    freepars=[FREE_PAR_1],
    initpoint=SYSTEM["ics"],
)
print("  Special points:", list_special_points(PC, "EP1"))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6 – CODIM-1: LIMIT-CYCLE BRANCH  (LC-C)
#             Attempts LC1 from every Hopf point found on EP1.
# ─────────────────────────────────────────────────────────────────────────────

for idx in range(1, 10):       # H1, H2, … (stops when point absent)
    sp_label = "H%d" % idx
    if _get_sp(PC, "EP1", sp_label) is None:
        break
    print("\n=== LC%d: limit cycle from EP1:%s ===" % (idx, sp_label))
    run_curve(
        PC, "LC%d" % idx, "LC-C",
        freepars=[FREE_PAR_1],
        initpoint="EP1:%s" % sp_label,
        overrides={"MaxNumPoints": 200, "LocBifPoints": ["PD", "LPC", "B"]},
    )
    print("  Special points on LC%d:" % idx, list_special_points(PC, "LC%d" % idx))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7 – CODIM-2: FOLD LOCUS  (LP-C)
#             One curve per LP point found on EP1.
#             Detects: Cusp (CP), Bogdanov–Takens (BT), Zero-Hopf (ZH).
# ─────────────────────────────────────────────────────────────────────────────

for idx in range(1, 10):
    sp_label = "LP%d" % idx
    if _get_sp(PC, "EP1", sp_label) is None:
        break
    print("\n=== LP_C%d: fold locus from EP1:%s ===" % (idx, sp_label))
    run_curve(
        PC, "LP_C%d" % idx, "LP-C",
        freepars=[FREE_PAR_1, FREE_PAR_2],
        initpoint="EP1:%s" % sp_label,
        overrides={
            "MaxNumPoints" : 150,
            "LocBifPoints" : ["CP", "BT", "ZH"],
        },
    )

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8 – CODIM-2: HOPF LOCUS  (H-C)
#             One curve per H point found on EP1.
#             Detects: BT, Generalized Hopf (GH), ZH, Double Hopf (HH).
# ─────────────────────────────────────────────────────────────────────────────

for idx in range(1, 10):
    sp_label = "H%d" % idx
    if _get_sp(PC, "EP1", sp_label) is None:
        break
    print("\n=== H_C%d: Hopf locus from EP1:%s ===" % (idx, sp_label))
    run_curve(
        PC, "H_C%d" % idx, "H-C",
        freepars=[FREE_PAR_1, FREE_PAR_2],
        initpoint="EP1:%s" % sp_label,
        overrides={
            "MaxNumPoints" : 150,
            "LocBifPoints" : ["BT", "GH", "ZH", "HH"],
        },
    )

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9 – PLOT
# ─────────────────────────────────────────────────────────────────────────────

PC.display(stability=True)
plt.title(SYSTEM["name"] + " — bifurcation diagram")
plt.xlabel(FREE_PAR_1)

savefig_path = os.getenv("PYDSTOOL_SAVEFIG")
if savefig_path:
    plt.savefig(savefig_path, dpi=160, bbox_inches="tight")
    print("Saved to", savefig_path)
else:
    plt.show()