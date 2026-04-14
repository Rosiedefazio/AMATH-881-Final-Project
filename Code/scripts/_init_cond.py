import numpy as np
from scipy.optimize import fsolve

# Define your system as a function that returns the RHS values
# The order of variables matters: C, T, M
def system_rhs(X, params):
    C, T, M = X
    r = params["r"]
    b = params["b"]
    gammaa = params["gammaa"]
    alphaa = params["alphaa"]
    lama = params["lama"]
    betaa = params["betaa"]
    deltaa = params["deltaa"]
    rhoa = params["rhoa"]
    etaa = params["etaa"]
    omegaa = params["omegaa"]

    dC_dt = r * C * (1 - b * C) * (C - gammaa) - alphaa * C * T
    dT_dt = lama * C + betaa * M * T - deltaa * T
    dM_dt = rhoa - etaa * M - omegaa * C * M
    return [dC_dt, dT_dt, dM_dt]

# Your provided parameters
pars = {
    "r": 0.08, "b": 0.01, "gammaa": 0, "alphaa": 0,
    "lama": 0.2, "betaa": 0.015900262162975, "deltaa": 1,
    "rhoa": 0, "etaa": 0, "omegaa": 1,
}
# Your provided initial guess
ics = {"C": 2, "T": 3, "M": 4}

# Convert initial conditions to a list/array in the same order as system_rhs expects
initial_guess_array = [ics["C"], ics["T"], ics["M"]]

# Use fsolve to find a root
# The 'args' parameter passes additional arguments to your system_rhs function
fixed_point_candidate = fsolve(system_rhs, initial_guess_array, args=(pars,))

print(f"Fixed point candidate: C={fixed_point_candidate[0]:.4f}, "
      f"T={fixed_point_candidate[1]:.4f}, M={fixed_point_candidate[2]:.4f}")

# Check if it's actually a root by plugging back into the RHS
rhs_at_candidate = system_rhs(fixed_point_candidate, pars)
print(f"RHS values at candidate: {rhs_at_candidate}")
if all(np.isclose(val, 0, atol=1e-6) for val in rhs_at_candidate):
    print("Candidate is a valid fixed point.")
else:
    print("Candidate is NOT a valid fixed point (RHS not close to zero).")