#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test script for pyDScont (Python 2.7)

Two example dynamical systems:
  * 3‑D Lorenz‑type model with 9 parameters
  * 2‑D predator‑prey model with 4 parameters

Each system is continued in one chosen parameter while the others are kept
fixed.  The script produces bifurcation diagrams using matplotlib.
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Import pyDScont – the exact import name may differ depending on how you
# installed it (e.g. `import pyDScont as ds` or `from pyDScont import ...`).
# Adjust the line below if necessary.
# ----------------------------------------------------------------------
try:
    import pyDScont as ds
except ImportError:
    raise ImportError(
        "pyDScont could not be imported.  Make sure the library is on your PYTHONPATH."
    )

# ----------------------------------------------------------------------
# Helper: generic continuation wrapper
# ----------------------------------------------------------------------
def run_continuation(system, init_state, param_names, param_values,
                     cont_param, cont_range, n_steps=200):
    """
    Perform a simple parameter continuation.

    Parameters
    ----------
    system          : callable(t, y, p) → dy/dt
    init_state      : list or np.ndarray, initial condition
    param_names     : list of strings, names of all parameters
    param_values    : list or np.ndarray, numerical values of all parameters
    cont_param      : string, name of the parameter to continue
    cont_range      : tuple (p_min, p_max)
    n_steps         : int, number of continuation points

    Returns
    -------
    param_vals, norm_vals : two 1‑D ndarrays
        * `param_vals` – the continuation parameter values
        * `norm_vals` – a scalar measure of the steady state (here ‖y‖)
    """
    # Build the pyDScont problem definition
    prob = ds.Problem()
    prob.set_ode(system, dim=len(init_state), param_names=param_names)

    # Fixed parameters (all except the continuation one)
    fixed_params = {name: val for name, val in zip(param_names, param_values)
                    if name != cont_param}
    prob.set_parameters(**fixed_params)

    # Initial guess for the equilibrium (use the supplied state)
    prob.set_initial_state(np.asarray(init_state, dtype=float))

    # Set up continuation
    cont = ds.Continuation(prob,
                          continuation_parameter=cont_param,
                          start=cont_range[0],
                          end=cont_range[1],
                          steps=n_steps,
                          method='pseudo_arclength')   # typical choice

    # Run the continuation
    cont.run()

    # Extract results (parameter value and Euclidean norm of the state)
    param_vals = np.array([pt[cont_param] for pt in cont.solution])
    norm_vals = np.linalg.norm(np.array([pt['state'] for pt in cont.solution]), axis=1)

    return param_vals, norm_vals


# ----------------------------------------------------------------------
# 1️⃣ 3‑D Lorenz‑type model (9 parameters)
# ----------------------------------------------------------------------
def lorenz_3d(t, y, p):
    """
    Generalized Lorenz system:

        dx/dt = sigma * (y - x) + a1 * x * z
        dy/dt = r * x - y - x * z + a2 * y * z
        dz/dt = -b * z + x * y + a3 * x**2 + a4 * y**2

    Parameters (9 total):
        sigma, r, b, a1, a2, a3, a4, c1, c2
    """
    x, y_, z = y
    sigma, r, b, a1, a2, a3, a4, c1, c2 = p

    dx = sigma * (y_ - x) + a1 * x * z + c1 * np.sin(z)
    dy = r * x - y_ - x * z + a2 * y_ * z + c2 * np.cos(x)
    dz = -b * z + x * y_ + a3 * x**2 + a4 * y_**2
    return np.array([dx, dy, dz], dtype=float)


# Parameter list for the 3‑D system
lorenz_params = [
    "sigma", "r", "b",          # classic Lorenz parameters
    "a1", "a2", "a3", "a4",    # nonlinear coupling coefficients
    "c1", "c2"                 # extra sinusoidal terms
]

# Choose a numeric set (feel free to modify)
lorenz_vals = [10.0, 28.0, 8.0/3.0,   # sigma, r, b
               0.1, 0.1, 0.05, 0.05,  # a1‑a4
               0.0, 0.0]              # c1, c2 (turned off for simplicity)

# Continuation in the Rayleigh number `r` from 0 to 60
lorenz_param_to_continue = "r"
lorenz_range = (0.0, 60.0)

# ----------------------------------------------------------------------
# 2️⃣ 2‑D Predator‑Prey model with harvesting (4 parameters)
# ----------------------------------------------------------------------
def predator_prey(t, y, p):
    """
    Modified Lotka‑Volterra equations with harvesting:

        dx/dt = a * x - b * x * y - h1 * x
        dy/dt = -c * y + d * x * y - h2 * y

    Parameters (4 total):
        a, b, c, d   – growth / interaction rates
        h1, h2       – harvesting rates (set to zero in the base case)
    """
    x, y_ = y
    a, b, c, d, h1, h2 = p
    dx = a * x - b * x * y_ - h1 * x
    dy = -c * y_ + d * x * y_ - h2 * y_
    return np.array([dx, dy], dtype=float)


# Parameter list for the 2‑D system (we keep the harvesting rates as extra entries)
pred_prey_params = ["a", "b", "c", "d", "h1", "h2"]

# Example numeric values (baseline Lotka‑Volterra, no harvesting)
pred_prey_vals = [1.0, 0.1, 1.5, 0.075,   # a, b, c, d
                  0.0, 0.0]               # h1, h2

# Continuation in the predator death rate `c` from 0.5 to 3.0
pred_prey_cont_param = "c"
pred_prey_range = (0.5, 3.0)

# ----------------------------------------------------------------------
# Main execution
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # --------------------------------------------------------------
    # 3‑D Lorenz continuation
    # --------------------------------------------------------------
    print("Running 3‑D Lorenz continuation (parameter = '{}')".format(lorenz_param_to_continue))
    lorenz_init = [1.0, 1.0, 1.0]   # a simple starting point
    lorenz_p, lorenz_norm = run_continuation(
        lorenz_3d,
        init_state=lorenz_init,
        param_names=lorenz_params,
        param_values=lorenz_vals,
        cont_param=lorenz_param_to_continue,
        cont_range=lorenz_range,
        n_steps=300
    )

    # Plot Lorenz bifurcation diagram
    plt.figure(figsize=(8, 5))
    plt.plot(lorenz_p, lorenz_norm, 'b.-')
    plt.title("Lorenz‑type system: norm of equilibrium vs. {}".format(lorenz_param_to_continue))
    plt.xlabel(lorenz_param_to_continue)
    plt.ylabel(r"$\| (x,y,z) \|_2$")
    plt.grid(True)

    # --------------------------------------------------------------
    # 2‑D Predator‑Prey continuation
    # --------------------------------------------------------------
    print("Running 2‑D predator‑prey continuation (parameter = '{}')".format(pred_prey_cont_param))
    pred_init = [10.0, 5.0]   # initial populations
    pred_p, pred_norm = run_continuation(
        predator_prey,
        init_state=pred_init,
        param_names=pred_prey_params,
        param_values=pred_prey_vals,
        cont_param=pred_prey_cont_param,
        cont_range=pred_prey_range,
        n_steps=250
    )

    # Plot predator‑prey bifurcation diagram
    plt.figure(figsize=(8, 5))
    plt.plot(pred_p, pred_norm, 'r.-')
    plt.title("Predator‑prey system: norm of equilibrium vs. {}".format(pred_prey_cont_param))
    plt.xlabel(pred_prey_cont_param)
    plt.ylabel(r"$\| (x,y) \|_2$")
    plt.grid(True)

    # --------------------------------------------------------------
    # Show both figures
    # --------------------------------------------------------------
    plt.show()
