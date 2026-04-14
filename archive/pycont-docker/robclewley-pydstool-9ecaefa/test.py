#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# --------------------------------------------------------------
# 1️⃣  Force a head‑less Matplotlib backend (no X‑server needed)
# --------------------------------------------------------------
import matplotlib
matplotlib.use('Agg')          # write to PNG files only

import matplotlib.pyplot as plt
import numpy as np

# --------------------------------------------------------------
# 2️⃣  PyDSTool imports
# --------------------------------------------------------------
try:
    from PyDSTool import args, Vode_ODEsystem, PyCont as pc
    from PyDSTool.errors import PyDSTool_KeyError
except ImportError as exc:
    raise ImportError(
        "PyDSTool could not be imported. Make sure it is installed and on PYTHONPATH."
    ) from exc


# --------------------------------------------------------------
# 3️⃣  Helper: generic continuation wrapper
# --------------------------------------------------------------
def run_continuation(
    varspecs,
    init_state,
    param_names,
    param_values,
    cont_param,
    cont_range,
    n_steps=200,
    model_name="model",
):
    """
    Perform a simple equilibrium continuation using PyDSTool 0.90.0.

    Returns
    -------
    param_vals : np.ndarray
        Continuation‑parameter values actually computed.
    norm_vals : np.ndarray
        Euclidean norm of the equilibrium point at each continuation step.
    """
    # ----- build DSargs -------------------------------------------------
    DSargs = args()
    DSargs.name = model_name
    DSargs.varspecs = varspecs
    DSargs.pars = dict(zip(param_names, param_values))

    var_names = list(varspecs.keys())
    DSargs.ics = {name: val for name, val in zip(var_names, init_state)}
    DSargs.tdomain = [0, 100]          # dummy – not used for equilibria

    # ----- ODE generator ------------------------------------------------
    ode = Vode_ODEsystem(DSargs)

    # ----- continuation set‑up (correct keyword names) ------------------
    cont_class = pc.ContClass(ode)

    curve_name = "EQ_curve"
    curve_initargs = {
        "name": curve_name,
        "type": "EP-C",               # equilibrium continuation
        "coordnames": var_names,      # state variables
        "parnames": [cont_param],     # parameter(s) to vary
        "freepars": [cont_param],     # same list – marked as “free”
    }
    cont_class.newCurve(curve_initargs)

    # ----- run the continuation -----------------------------------------
    step = (cont_range[1] - cont_range[0]) / float(n_steps)

    curve = cont_class[curve_name]
    curve.StepSize = step
    curve.maxNumPoints = n_steps
    curve.forward()

    # ----- extract results -----------------------------------------------
    sol = curve.sol

    # continuation‑parameter values (already a 1‑D NumPy array)
    param_vals = np.asarray(sol[cont_param], dtype=float)

    # equilibrium points – try the old “points” entry first
    try:
        points = sol["points"].coordarray.T          # (n_points, n_vars)
    except PyDSTool_KeyError:
        points = np.column_stack(
            [np.asarray(sol[v], dtype=float) for v in var_names]
        )

    norm_vals = np.linalg.norm(points, axis=1)

    return param_vals, norm_vals


# ----------------------------------------------------------------------
# 1️⃣ 3‑D Lorenz‑type model (9 parameters)
# ----------------------------------------------------------------------
lorenz_varspecs = {
    "x": "sigma*(y - x) + a1*x*z + c1*sin(z)",
    "y": "r*x - y - x*z + a2*y*z + c2*cos(x)",
    "z": "-b*z + x*y + a3*x**2 + a4*y**2",
}
lorenz_params = [
    "sigma", "r", "b", "a1", "a2", "a3", "a4", "c1", "c2",
]
lorenz_vals = [
    10.0, 28.0, 8.0 / 3.0, 0.1, 0.1, 0.05, 0.05, 0.0, 0.0,
]
lorenz_cont_param = "r"
lorenz_range = (0.0, 60.0)


# ----------------------------------------------------------------------
# 2️⃣ 2‑D Predator–Prey model with harvesting (6 parameters)
# ----------------------------------------------------------------------
predator_varspecs = {
    "x": "a*x - b*x*y - h1*x",
    "y": "-c*y + d*x*y - h2*y",
}
predator_params = ["a", "b", "c", "d", "h1", "h2"]
predator_vals = [1.0, 0.1, 1.5, 0.075, 0.0, 0.0]
predator_cont_param = "c"
predator_range = (0.5, 3.0)


# ----------------------------------------------------------------------
# Main execution
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # -------------------- Lorenz continuation --------------------
    print(
        "Running 3‑D Lorenz continuation (parameter = '{}')".format(
            lorenz_cont_param
        )
    )
    lorenz_init = [1.0, 1.0, 1.0]

    lorenz_p, lorenz_norm = run_continuation(
        lorenz_varspecs,
        init_state=lorenz_init,
        param_names=lorenz_params,
        param_values=lorenz_vals,
        cont_param=lorenz_cont_param,
        cont_range=lorenz_range,
        n_steps=300,
        model_name="lorenz",
    )

    plt.figure(figsize=(8, 5))
    plt.plot(lorenz_p, lorenz_norm, "b.-")
    plt.title(
        "Lorenz‑type system: norm of equilibrium vs. {}".format(lorenz_cont_param)
    )
    plt.xlabel(lorenz_cont_param)
    plt.ylabel(r"$\|(x,y,z)\|_2$")
    plt.grid(True)
    plt.savefig("lorenz_norm.png")
    plt.close()

    # -------------------- Predator‑prey continuation --------------------
    print(
        "Running 2‑D predator‑prey continuation (parameter = '{}')".format(
            predator_cont_param
        )
    )
    pred_init = [10.0, 5.0]

    pred_p, pred_norm = run_continuation(
        predator_varspecs,
        init_state=pred_init,
        param_names=predator_params,
        param_values=predator_vals,
        cont_param=predator_cont_param,
        cont_range=predator_range,
        n_steps=250,
        model_name="predator",
    )

    plt.figure(figsize=(8, 5))
    plt.plot(pred_p, pred_norm, "r.-")
    plt.title(
        "Predator‑prey system: norm of equilibrium vs. {}".format(
            predator_cont_param
        )
    )
    plt.xlabel(predator_cont_param)
    plt.ylabel(r"$\|(x,y)\|_2$")
    plt.grid(True)
    plt.savefig("predator_norm.png")
    plt.close()

    print("Figures saved as 'lorenz_norm.png' and 'predator_norm.png'")