# 3-DoF Longitudinal Aircraft System Identification

A MATLAB implementation of a classical system identification pipeline for estimating longitudinal aerodynamic parameters of a fixed-wing aircraft from flight data. This serves as a **baseline reference implementation** — the results it produces are intended to provide a basis of comparison for other system identification methods.

The workflow uses two sequential stages: an **Equation-Error** method for initial parameter estimates, followed by an **Output-Error** method for maximum-likelihood refinement.

---

## Files

| File | Description |
|---|---|
| `AC_VER2_main.m` | Script to call the other files |
| `AC_VER2_sysID.m` | Main sysID implementation — runs both equation-error and output-error methods |
| `AC_VER2_statespace.m` | Supporting file defining the nonlinear 3-DoF fixed-wing state-space model |
| `AC_VER2_statespace.m` | Supporting file for plotting states and inputs |
| `model.py` / `classical_OEM_resultsforAIAA 1` | Files used to generate the flight data |

---

## Pipeline

```
Flight Data
    │
    ▼
Equation-Error Method  —  constrained least squares (lsqlin)
    │                     decomposes into lift, drag, and moment sub-problems
    ▼  initial parameter estimates
Output-Error Method  —  oe.m (SIDPAC)
    │                   Modified Newton-Raphson + Simplex fallback
    │                   Full nonlinear RK4 simulation at each iteration
    ▼
Final Parameter Estimates + Cramér-Rao Bounds
```

---

## Parameters Identified

```
CL0, CL_alpha               — lift model
CD0, k_CD                   — drag polar
CM0, CM_alpha, CM_Q, CM_e   — pitching moment model
```

---

## Dependencies

- MATLAB + Optimization Toolbox (`lsqlin`)
- [SIDPAC](https://software.nasa.gov/) toolbox (`oe.m` and supporting routines: `mnr.m`, `simplex.m`, `estrr.m`, `misvd.m`, `cvec.m`, `compcost.m`)
- `deriv.m` — smoothed numerical differentiation
- `lesq.m` — unconstrained least squares
