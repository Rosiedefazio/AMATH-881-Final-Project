"""
Bifurcation analysis of a cancer–immune interaction model using PyDSTool/PyCont.

Codim-1 curves:
  EP-C  — equilibrium point continuation (detects LP/fold and H/Hopf points)
  LC-C  — limit cycle continuation (from Hopf point)

Codim-2 curves:
  LP-C  — fold (saddle-node) locus in 2-parameter plane
  H-C2  — Hopf locus in 2-parameter plane
"""

import PyDSTool as dst
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ============================================================
# 0. VERSION-AGNOSTIC HELPERS
#    Attribute-only — no getSpecialPoints() call.
#    Tested against the classic PyDSTool/PyCont storage layout
#    where special points live in a dict keyed by their
#    registration name ("H1", "LP1", etc.).
# ============================================================
# ============================================================
# 0. VERSION-AGNOSTIC HELPERS
#    Attribute-only — no getSpecialPoints() call.
#
#    KEY INSIGHT: ContClass.newCurve validates initpoint names
#    against curve.BifPoints, so BifPoints is checked FIRST.
#    specialpoints may use a different (shorter) labelling
#    scheme and must not take priority.
# ============================================================
# ============================================================
# 0. HELPERS
# ============================================================

def _point_label(pt):
    for attr in ("label", "name"):
        val = getattr(pt, attr, None)
        if val is not None:
            return str(val)
    return repr(pt)


def dump_curve_internals(curve, curve_name):
    """
    Print every attribute on the curve that could hold special points.
    Run this once after a forward()/backward() call to see the real
    keys your PyDSTool build uses.
    """
    print("\n==== RAW INTERNALS: '{}' ====".format(curve_name))
    for attr in ("BifPoints", "specialpoints", "curves", "labels"):
        val = getattr(curve, attr, "** ATTRIBUTE MISSING **")
        print("  .{:20s} -> {!r}".format(attr, val))

    # Also walk the object dict for anything else that looks relevant
    for k, v in vars(curve).items():
        if k not in ("BifPoints", "specialpoints", "curves", "labels"):
            if "point" in k.lower() or "bif" in k.lower():
                print("  .{:20s} -> {!r}".format(k, v))
    print("=" * 40)


def get_valid_initpoint_keys(curve, curve_name):
    """
    Return the dict that ContClass.newCurve validates initpoint names
    against, along with the source attribute name found.

    ContClass.newCurve does (roughly):
        pointname = initpoint.split(':')[1]
        if pointname not in self[curvename].BifPoints: raise KeyError

    So we must return keys from the exact same dict it uses.
    """
    # Try every plausible attribute name in the order ContClass checks them
    for attr in ("BifPoints", "specialpoints", "LocBifPoints"):
        sp = getattr(curve, attr, None)
        if isinstance(sp, dict) and sp:
            print("  [get_valid_initpoint_keys] using '{}': keys = {}".format(
                attr, list(sp.keys())
            ))
            return sp

    # Nothing dict-like found — return empty
    print("  [get_valid_initpoint_keys] no dict attribute found on curve")
    return {}


def print_special_points(PC, curve_name):
    curve = PC[curve_name]
    kv    = get_valid_initpoint_keys(curve, curve_name)

    print("\n--- Special points on curve '{}' ---".format(curve_name))
    if not kv:
        print("  (none detected)")
        return

    for key, pt in kv.items():
        try:
            vals = {k: "{}".format(v) for k, v in pt.coorddict.items()}
        except AttributeError:
            vals = "(coordinate dict unavailable)"
        print("  [{}]  {}".format(key, vals))
# ============================================================
# 1. MODEL DEFINITION
# ============================================================

pars = {
    "r"     : 0.50,
    "b"     : 0.002,
    "gammaa": 0.10,
    "alphaa": 0.50,
    "lama"  : 0.01,
    "betaa" : 0.10,
    "deltaa": 0.20,
    "rhoa"  : 0.10,
    "etaa"  : 0.05,
    "omegaa": 0.01,
}

ics = {
    "C": 0.50,
    "T": 0.30,
    "M": 1.50,
}

varspecs = {
    "C": "r*C*(1 - b*C)*(C - gammaa) - alphaa*C*T",
    "T": "lama*C + betaa*M*T - deltaa*T",
    "M": "rhoa - etaa*M - omegaa*C*M",
}

myDSargs          = dst.args(name="cancer_immune")
myDSargs.pars     = pars
myDSargs.varspecs = varspecs
myDSargs.ics      = ics
myDSargs.tdomain  = [0, 500]

myode = dst.Generator.Vode_ODEsystem(myDSargs)

# ============================================================
# 2. CONTINUATION OBJECT
# ============================================================

PC = dst.ContClass(myode)

# ============================================================
# 3. CODIM-1 — Equilibrium continuation in alpha
# ============================================================

print("=" * 60)
print("Codim-1: EP-C  — alpha as free parameter")
print("=" * 60)

PCargs              = dst.args(name="EQ_alpha", type="EP-C")
PCargs.freepars     = ["alphaa"]
PCargs.MaxNumPoints = 400
PCargs.MaxStepSize  = 0.02
PCargs.MinStepSize  = 1e-6
PCargs.StepSize     = 5e-3
PCargs.LocBifPoints = "all"
PCargs.SaveEigen    = True
PCargs.verbosity    = 2

PC.newCurve(PCargs)
PC["EQ_alpha"].forward()
PC["EQ_alpha"].backward()

print_special_points(PC, "EQ_alpha")

# ============================================================
# 4. CODIM-1 — Equilibrium continuation in rho
# ============================================================

print("\n" + "=" * 60)
print("Codim-1: EP-C  — rho as free parameter")
print("=" * 60)

PCargs              = dst.args(name="EQ_rho", type="EP-C")
PCargs.freepars     = ["rhoa"]
PCargs.MaxNumPoints = 400
PCargs.MaxStepSize  = 0.02
PCargs.MinStepSize  = 1e-6
PCargs.StepSize     = 5e-3
PCargs.LocBifPoints = "all"
PCargs.SaveEigen    = True
PCargs.verbosity    = 2

PC.newCurve(PCargs)
PC["EQ_rho"].forward()
PC["EQ_rho"].backward()

print_special_points(PC, "EQ_rho")

# ============================================================
# 5. CODIM-1 — Limit cycle continuation from Hopf on EQ_alpha
# ============================================================

# Dump raw internals ONCE so you can see the real registered key names
dump_curve_internals(PC["EQ_alpha"], "EQ_alpha")

# Pull the authoritative dict ContClass.newCurve validates against
eq_alpha_bif = get_valid_initpoint_keys(PC["EQ_alpha"], "EQ_alpha")

# Hopf keys: start with H but not GH or ZH
hopf_keys = [
    k for k in eq_alpha_bif
    if k.upper().startswith("H") and not k.upper().startswith(("GH", "ZH"))
]

# Fold keys
fold_keys = [
    k for k in eq_alpha_bif
    if k.upper().startswith("LP")
]

print("\nHopf keys found : {}".format(hopf_keys))
print("Fold keys found : {}".format(fold_keys))

hopf_key = hopf_keys[0] if hopf_keys else None
fold_key  = fold_keys[0] if fold_keys else None

if hopf_key:
    print("\n" + "=" * 60)
    print("Codim-1: LC-C  — from {} on EQ_alpha".format(hopf_key))
    print("=" * 60)

    PCargs              = dst.args(name="LC_alpha", type="LC-C")
    PCargs.initpoint    = "EQ_alpha:{}".format(hopf_key)
    PCargs.freepars     = ["alphaa"]
    PCargs.MaxNumPoints = 300
    PCargs.MaxStepSize  = 0.05
    PCargs.MinStepSize  = 1e-6
    PCargs.StepSize     = 0.01
    PCargs.LocBifPoints = ["PD", "LPC"]
    PCargs.SaveEigen    = True
    PCargs.verbosity    = 2

    PC.newCurve(PCargs)
    PC["LC_alpha"].forward()
    PC["LC_alpha"].backward()
    print_special_points(PC, "LC_alpha")

else:
    print(
        "\n[!] No Hopf key found in BifPoints — skipping LC continuation.\n"
        "    Check the dump above for the real registered key names."
    )

# ============================================================
# 6. CODIM-2 — Hopf locus in (alpha, rho) plane (H-C2)
# ============================================================

if hopf_key:
    print("\n" + "=" * 60)
    print("Codim-2: H-C2  — Hopf curve in (alpha, rho)")
    print("=" * 60)

    PCargs              = dst.args(name="Hopf2D", type="H-C2")
    PCargs.initpoint    = "EQ_alpha:{}".format(hopf_key)
    PCargs.freepars     = ["alphaa", "rhoa"]
    PCargs.MaxNumPoints = 300
    PCargs.MaxStepSize  = 0.05
    PCargs.MinStepSize  = 1e-6
    PCargs.StepSize     = 0.01
    PCargs.LocBifPoints = ["GH", "BT", "ZH"]
    PCargs.SaveEigen    = True
    PCargs.verbosity    = 2

    PC.newCurve(PCargs)
    PC["Hopf2D"].forward()
    PC["Hopf2D"].backward()
    print_special_points(PC, "Hopf2D")

# ============================================================
# 7. CODIM-2 — Fold (LP) locus in (alpha, rho) plane (LP-C)
# ============================================================

if fold_key:
    print("\n" + "=" * 60)
    print("Codim-2: LP-C  — Fold curve in (alpha, rho)")
    print("=" * 60)

    PCargs              = dst.args(name="Fold2D", type="LP-C")
    PCargs.initpoint    = "EQ_alpha:{}".format(fold_key)
    PCargs.freepars     = ["alphaa", "rhoa"]
    PCargs.MaxNumPoints = 300
    PCargs.MaxStepSize  = 0.05
    PCargs.MinStepSize  = 1e-6
    PCargs.StepSize     = 0.01
    PCargs.LocBifPoints = ["CP", "BT", "ZH"]
    PCargs.SaveEigen    = True
    PCargs.verbosity    = 2

    PC.newCurve(PCargs)
    PC["Fold2D"].forward()
    PC["Fold2D"].backward()
    print_special_points(PC, "Fold2D")

else:
    print(
        "\n[!] No fold key found in BifPoints — skipping LP-C curve.\n"
        "    Check the dump above for the real registered key names."
    )

# ============================================================
# 8. PLOTTING
# ============================================================

fig = plt.figure(figsize=(16, 12))
gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.35)

# --- 8a: Codim-1 alpha vs C --------------------------------------------
ax1 = fig.add_subplot(gs[0, 0])
PC["EQ_alpha"].display(["alphaa", "C"], stability=True, axes=ax1)
if "LC_alpha" in PC.curves:
    PC["LC_alpha"].display(["alphaa", "C"], stability=True, axes=ax1)
ax1.set_title("Codim-1: alpha  ->  C\n(equilibria + limit cycles)")
ax1.set_xlabel("alpha  (kill rate)")
ax1.set_ylabel("C  (cancer cells)")
ax1.grid(True, linestyle="--", alpha=0.5)
ax1.legend(
    ["Stable eq.", "Unstable eq.", "Stable LC", "Unstable LC"],
    fontsize=8,
)

# --- 8b: Codim-1 rho vs C ----------------------------------------------
ax2 = fig.add_subplot(gs[0, 1])
PC["EQ_rho"].display(["rhoa", "C"], stability=True, axes=ax2)
ax2.set_title("Codim-1: rho  ->  C\n(equilibria)")
ax2.set_xlabel("rho  (M-cell influx)")
ax2.set_ylabel("C  (cancer cells)")
ax2.grid(True, linestyle="--", alpha=0.5)

# --- 8c: Codim-2 Hopf locus --------------------------------------------
ax3 = fig.add_subplot(gs[1, 0])
if "Hopf2D" in PC.curves:
    PC["Hopf2D"].display(["alphaa", "rhoa"], axes=ax3)
    ax3.set_title("Codim-2: Hopf curve\nin (alpha, rho) plane")
    ax3.set_xlabel("alpha")
    ax3.set_ylabel("rho")
    ax3.grid(True, linestyle="--", alpha=0.5)
else:
    ax3.text(
        0.5, 0.5,
        "No Hopf curve\n(no H point detected)",
        ha="center", va="center", transform=ax3.transAxes,
    )
    ax3.set_title("Codim-2: Hopf curve (unavailable)")

# --- 8d: Codim-2 Fold locus --------------------------------------------
ax4 = fig.add_subplot(gs[1, 1])
if "Fold2D" in PC.curves:
    PC["Fold2D"].display(["alphaa", "rhoa"], axes=ax4)
    ax4.set_title("Codim-2: Fold (saddle-node) curve\nin (alpha, rho) plane")
    ax4.set_xlabel("alpha")
    ax4.set_ylabel("rho")
    ax4.grid(True, linestyle="--", alpha=0.5)
else:
    ax4.text(
        0.5, 0.5,
        "No fold curve\n(no LP point detected)",
        ha="center", va="center", transform=ax4.transAxes,
    )
    ax4.set_title("Codim-2: Fold curve (unavailable)")

fig.suptitle(
    "Cancer-Immune ODE  —  PyCont Bifurcation Diagrams",
    fontsize=14,
    fontweight="bold",
)
plt.savefig("bifurcation_diagrams.png", dpi=150, bbox_inches="tight")
print("\nDone. Diagrams saved to bifurcation_diagrams.png")