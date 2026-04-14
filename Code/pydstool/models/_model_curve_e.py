"""
Template: Codim-1 and Codim-2 bifurcation analysis for 3D ODEs
using PyDSTool (PyDSCont).

Workflow (queue-driven, recursive branch following)
----------------------------------------------------
  EP-C  →  LP → LP-C (codim-2 fold locus)
        →  H  → LC-C  (limit-cycle branch)
               → H-C  (codim-2 Hopf locus)
        →  BP → new EP-C  (secondary equilibrium branch)
  LC-C  →  PD  → new LC-C  (period-doubled branch)
         →  LPC → LP-C     (codim-2 fold-of-cycles locus)

The BFS queue ensures every discovered special point seeds
exactly one new curve (duplicates are suppressed).
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from collections import deque
from time import perf_counter
from typing import Dict, List, Optional, Any

from PyDSTool import *

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1 – SYSTEM DEFINITION
# ─────────────────────────────────────────────────────────────────────────────

SYSTEM = dict(
    name="BiologicalModel",

    pars={
        "r": 0.08,
        "b": 0.01,
        "gammaa": 0,
        "alphaa": 0,
        "lama": 0.2,
        "betaa": 0.015900262162975,
        "deltaa": 1,
        "rhoa": 0,
        "etaa": 0,
        "omegaa": 1,
    },

    varspecs={
        "C": "r*C*(1 - b*C)*(C - gammaa) - alphaa*C*T",
        "T": "lama*C + betaa*M*T - deltaa*T",
        "M": "rhoa - etaa*M - omegaa*C*M",
    },

    ics={"C": 8.3, "T": 2.3, "M": 6.9},

    tdata=[-50000000.0, 50000000000000.0],
)

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2 – CONTINUATION SETTINGS
# ─────────────────────────────────────────────────────────────────────────────

FREE_PAR_1 = "alphaa"
FREE_PAR_2 = "gammaa"

PAR_DOMAIN = {
    FREE_PAR_1: [0.0, 0.5],
    FREE_PAR_2: [0.0, 30],
}

BASE_CONT = dict(
    StepSize     = 1e-2,
    MaxStepSize  = 5e-1,
    MaxNumPoints = 300,
    verbosity    = 2,
    SaveEigen    = True,
    LocBifPoints = "all",
    StopAtPoints = ["B"],
)

# Safety cap – prevents runaway cascades on highly branched systems.
MAX_CURVES = 60

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3 – HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def build_ode(sys_dict: dict) -> Any:
    """Construct a PyDSTool Vode_ODEsystem from a plain dict."""
    d = args(name=sys_dict["name"])
    d.pars     = dict(sys_dict["pars"])
    d.varspecs = dict(sys_dict["varspecs"])
    d.ics      = dict(sys_dict["ics"])
    d.tdata    = list(sys_dict["tdata"])
    d.pdomain  = {k: [-1e9, 1e9] for k in sys_dict["pars"]}
    return Generator.Vode_ODEsystem(d)


def _get_sp(PC, curve: str, label: str) -> Any:
    try:
        return PC[curve].getSpecialPoint(label)
    except Exception:
        return None


def list_special_points(PC, curve: str) -> List[str]:
    try:
        return sorted(PC[curve].getSpecialPoint().keys())
    except Exception:
        return []


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
    """Register, run (fwd + optional bwd), and clean a continuation curve."""
    s = dict(BASE_CONT)
    if overrides:
        s.update(overrides)

    a = args(name=name, type=curve_type)
    a.freepars  = freepars
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


# ── Name counters & duplicate-seed guard ─────────────────────────────────────

_name_counters: Dict[str, int] = {}
_used_seeds: set = set()


def _next_name(prefix: str) -> str:
    """Return the next unique curve name for the given prefix."""
    _name_counters[prefix] = _name_counters.get(prefix, 0) + 1
    return "%s%d" % (prefix, _name_counters[prefix])


def _claim_seed(source_curve: str, sp_label: str) -> bool:
    """
    Mark (source_curve, sp_label) as consumed.
    Returns True the first time, False on every subsequent call.
    Prevents the same special point from spawning duplicate curves.
    """
    key = "%s:%s" % (source_curve, sp_label)
    if key in _used_seeds:
        return False
    _used_seeds.add(key)
    return True


# ── Branch-discovery logic ───────────────────────────────────────────────────

def discover_branches(
    PC,
    curve_name: str,
    curve_type: str,
    codim2: bool = True,
) -> List[Dict]:
    """
    Inspect special points on a completed curve and return a list of
    new work-queue task dicts.

    Recognised bifurcations per curve type
    ───────────────────────────────────────
    EP-C  →  H   : LC-C  (codim-1 limit cycle)
                   H-C   (codim-2 Hopf locus)         [if codim2]
          →  LP  : LP-C  (codim-2 fold locus)         [if codim2]
          →  BP  : EP-C  (secondary equilibrium branch)

    LC-C  →  PD  : LC-C  (period-doubled branch)
          →  LPC : LP-C  (codim-2 fold-of-cycles)     [if codim2]
    """
    sps = list_special_points(PC, curve_name)
    print("  Special points on %s: %s" % (curve_name, sps if sps else "none"))

    tasks: List[Dict] = []

    for sp in sps:
        btype = sp.rstrip("0123456789")   # strip trailing digits → "H","LP","BP"…

        # ── Equilibrium branch ──────────────────────────────────────────────
        if curve_type == "EP-C":

            if btype == "H":
                if not _claim_seed(curve_name, sp):
                    continue
                # Codim-1: limit-cycle branch
                tasks.append(dict(
                    name      = _next_name("LC"),
                    type      = "LC-C",
                    freepars  = [FREE_PAR_1],
                    initpoint = "%s:%s" % (curve_name, sp),
                    overrides = {
                        "MaxNumPoints": 200,
                        "LocBifPoints": ["PD", "LPC", "B"],
                    },
                ))
                # Codim-2: Hopf locus
                if codim2:
                    tasks.append(dict(
                        name      = _next_name("H_C"),
                        type      = "H-C",
                        freepars  = [FREE_PAR_1, FREE_PAR_2],
                        initpoint = "%s:%s" % (curve_name, sp),
                        overrides = {
                            "MaxNumPoints": 150,
                            "LocBifPoints": ["BT", "GH", "ZH", "HH"],
                        },
                    ))

            elif btype == "LP":
                if not _claim_seed(curve_name, sp):
                    continue
                # Codim-2: fold locus
                if codim2:
                    tasks.append(dict(
                        name      = _next_name("LP_C"),
                        type      = "LP-C",
                        freepars  = [FREE_PAR_1, FREE_PAR_2],
                        initpoint = "%s:%s" % (curve_name, sp),
                        overrides = {
                            "MaxNumPoints": 150,
                            "LocBifPoints": ["CP", "BT", "ZH"],
                        },
                    ))

            elif btype == "BP":
                if not _claim_seed(curve_name, sp):
                    continue
                # Codim-1: secondary equilibrium branch
                tasks.append(dict(
                    name      = _next_name("EP"),
                    type      = "EP-C",
                    freepars  = [FREE_PAR_1],
                    initpoint = "%s:%s" % (curve_name, sp),
                    overrides = {},
                ))

        # ── Limit-cycle branch ──────────────────────────────────────────────
        elif curve_type == "LC-C":

            if btype == "PD":
                if not _claim_seed(curve_name, sp):
                    continue
                # Codim-1: period-doubled branch
                tasks.append(dict(
                    name      = _next_name("LC"),
                    type      = "LC-C",
                    freepars  = [FREE_PAR_1],
                    initpoint = "%s:%s" % (curve_name, sp),
                    overrides = {
                        "MaxNumPoints": 200,
                        "LocBifPoints": ["PD", "LPC", "B"],
                    },
                ))

            elif btype == "LPC":
                if not _claim_seed(curve_name, sp):
                    continue
                # Codim-2: fold-of-cycles locus
                if codim2:
                    tasks.append(dict(
                        name      = _next_name("LPC_C"),
                        type      = "LP-C",
                        freepars  = [FREE_PAR_1, FREE_PAR_2],
                        initpoint = "%s:%s" % (curve_name, sp),
                        overrides = {
                            "MaxNumPoints": 150,
                            "LocBifPoints": ["CP", "BT", "ZH"],
                        },
                    ))

        # Codim-2 curves (H-C, LP-C) are terminal: their own special points
        # (BT, GH, ZH, HH, CP) are co-dimension-2 and cannot seed further
        # standard continuation curves in PyDSTool.

    return tasks


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4 – BUILD SYSTEM
# ─────────────────────────────────────────────────────────────────────────────

ode = build_ode(SYSTEM)
PC  = ContClass(ode)

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5 – BFS QUEUE-DRIVEN CONTINUATION
# ─────────────────────────────────────────────────────────────────────────────

# Seed the queue with the primary equilibrium branch.
work_queue: deque = deque([
    dict(
        name      = _next_name("EP"),   # → "EP1"
        type      = "EP-C",
        freepars  = [FREE_PAR_1],
        initpoint = SYSTEM["ics"],
        overrides = {},
    )
])

completed:  List[str] = []
skipped:    List[str] = []
curve_count = 0

print("\n" + "=" * 70)
print("  Starting BFS bifurcation continuation (cap = %d curves)" % MAX_CURVES)
print("=" * 70)

while work_queue and curve_count < MAX_CURVES:
    task = work_queue.popleft()
    curve_count += 1

    name      = task["name"]
    ctype     = task["type"]
    freepars  = task["freepars"]
    initpoint = task["initpoint"]
    overrides = task.get("overrides", {})

    print(
        "\n[%d/%d] %s (%s)  initpoint=%s"
        % (curve_count, MAX_CURVES, name, ctype,
           initpoint if isinstance(initpoint, str) else "ICs")
    )

    ok = run_curve(
        PC, name, ctype,
        freepars  = freepars,
        initpoint = initpoint,
        overrides = overrides,
    )

    if ok:
        completed.append(name)
        # Discover and enqueue all new branches.
        new_tasks = discover_branches(PC, name, ctype, codim2=True)
        if new_tasks:
            print(
                "  → queuing %d new curve(s): %s"
                % (len(new_tasks), [t["name"] for t in new_tasks])
            )
        work_queue.extend(new_tasks)
    else:
        skipped.append(name)

if work_queue:
    remaining = [t["name"] for t in work_queue]
    print(
        "\n⚠  Safety cap (%d) reached. Unrun curves: %s"
        % (MAX_CURVES, remaining)
    )

print("\n" + "=" * 70)
print("  Completed : %s" % completed)
print("  Skipped   : %s" % skipped)
print("=" * 70)

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6 – PLOT
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